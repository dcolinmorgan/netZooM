chipdir='ScanBedResults/'; % directory with input bed files
outdir='PromoMapRefSeq_100kb/'; % directory to print mapped files (make sure this is already created)
chipfiles=dir([chipdir, '*.bed']);

tssmaxdist=100000; % region around tss to print
pvalcutoff=1e-4; % only keep motif hits with this significance
% parameters for processing a subset of the bed files in chipdir
mindix=1; 
maxidx=10;
% maxidx=length(chipfiles); % to run all bed files

% read in the annotation information
fid=fopen('./ReferenceData/refseq_hg19_11162016', 'r');
RefGene=textscan(fid, '%u%s%s%s%f%f%f%f%f%s%s%s%s%s%s%s', 'delimiter', '\t', 'headerlines', 1);
fclose(fid);
% only keep genes on the canonical chromosomes
allchr={'chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr20','chr21','chr22','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrM','chrX','chrY'};
filter=ismember(RefGene{3}, allchr);
for(cnt=1:length(RefGene))
	RefGene{cnt}=RefGene{cnt}(filter);
end

% break up annotation information by chromosome
mcnt=1;
uchr=unique(RefGene{3});
chridx=cell(length(uchr),1);
chrTSS=cell(length(uchr),1);
chrStrand=cell(length(uchr),1);
chrStrandBool=cell(length(uchr),1);
for(cnt=1:length(uchr))
        chridx{cnt}=find(strcmp(RefGene{3}, uchr{cnt}));
        locStrand=RefGene{4}(chridx{cnt});
        locEnd=RefGene{6}(chridx{cnt});
        chrTSS{cnt}=RefGene{5}(chridx{cnt});
        chrTSS{cnt}(strcmp(locStrand, '-'))=locEnd(strcmp(locStrand, '-'));
        chrTSS{cnt}=double(chrTSS{cnt});
        chrStrand{cnt}=locStrand;
        chrStrandBool{cnt}=ones(length(locStrand),1);
	chrStrandBool{cnt}(strcmp(locStrand, '-'))=-1*chrStrandBool{cnt}(strcmp(locStrand, '-'));

	[uStart,sidx,sloc]=unique(chrTSS{cnt});
	keepidx=zeros(length(sidx),1);

	locGene=RefGene{13}(chridx{cnt});
	for(scnt=1:length(uStart))
		if(length(unique(locGene(sloc==scnt)))==1)
			keepidx(scnt)=1;
		else
			mcnt=mcnt+1;
		end
	end
	chridx{cnt}=chridx{cnt}(sidx(keepidx==1));
	chrTSS{cnt}=chrTSS{cnt}(sidx(keepidx==1));
	chrStrand{cnt}=chrStrand{cnt}(sidx(keepidx==1));
	chrStrandBool{cnt}=chrStrandBool{cnt}(sidx(keepidx==1));
end
RefGene=RefGene{13};

% iterate through TF-motif bed files and map to genes
for(cnt=mindix:1:maxidx)
	tic
	fid=fopen([chipdir, chipfiles(cnt).name], 'r');
	A=textscan(fid, '%s%f%f%s%f%s%f%f%f%f', 'delimiter', '\t');
	fclose(fid);
	filter=ismember(A{1}, uchr) & A{7}<=pvalcutoff; % filter which peaks to map
	if(sum(~filter))
		for(c=1:length(A))
			A{c}=A{c}(filter);
		end
	end
	peaksummit=(A{2}+A{3})/2;
	TFName=chipfiles(cnt).name(1:end-4);
	disp(['Now working on ', TFName, '......']);
	[~,locchr]=ismember(A{1}, uchr);
	LocGene=cell(length(locchr),1);
	DistVals=zeros(length(locchr),1);
	if(chipfiles(cnt).bytes<1.5e9)
		% iterate by chromosome (maps all the peaks on the same chromosome at the same time)
		for(chrcnt=1:length(uchr))
			temp1=repmat(chrTSS{chrcnt}, 1, sum(locchr==chrcnt));
			temp2=repmat(peaksummit(locchr==chrcnt)', length(chrTSS{chrcnt}),1);
			[~,idx]=min(abs(temp1-temp2));
			DistVals(locchr==chrcnt)=chrStrandBool{chrcnt}(idx).*(peaksummit(locchr==chrcnt)-chrTSS{chrcnt}(idx));
			LocGene(locchr==chrcnt)=RefGene(chridx{chrcnt}(idx));
		end
	else
		% iterate by peak (maps peaks one-by-one)
		for(peakcnt=1:length(peaksummit))
			[~,idx]=min(abs(chrTSS{locchr(peakcnt)}-peaksummit(peakcnt)));
			DistVals(peakcnt)=chrStrandBool{locchr(peakcnt)}(idx)*(peaksummit(peakcnt)-chrTSS{locchr(peakcnt)}(idx));
			LocGene{peakcnt}=RefGene{chridx{locchr(peakcnt)}(idx)};
		end
	end
	filter=abs(DistVals)<tssmaxdist; % filter which peaks to print
	fid=fopen([outdir, TFName,'_map.txt'], 'W');
	for(peakcnt=1:length(peaksummit))
		if(filter(peakcnt))
			fprintf(fid, '%s\t%.4f\t%g\t%s\t%.1f\n', [A{1}{peakcnt}, ':', num2str(A{2}(peakcnt)), '-', num2str(A{3}(peakcnt))], A{5}(peakcnt), A{7}(peakcnt), LocGene{peakcnt}, DistVals(peakcnt));
		end
	end
	fclose(fid);
	timelapse=toc;
	disp(['This motif took ', num2str(timelapse), ' seconds to map.']);
end
