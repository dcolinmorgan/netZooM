PWMdir='../tfprior/cisbp/pwms_all_motifs/';
outdir='MotifsInUniprobe/';
TFlist=textread('../tfprior/cisbp/directandinferred_uniqueTF_infcon.txt', '%s', 'delimiter', '\n');
% filelist=dir([PWMdir, '*.txt']);

for(cnt=1:length(TFlist))
% for(cnt=1:length(filelist))
	% TF=filelist(cnt).name(1:end-4);
	TF=TFlist{cnt};
	disp(['Now working on ', TF, '......']);
	[Pos, A, C, G, T]=textread([PWMdir, TF, '.txt'], '%u%f%f%f%f', 'headerlines', 1);
	if(length(Pos))
		fid=fopen([outdir, TF, '.txt'], 'wt');
		fprintf(fid, '%s\n', TF);
		fprintf(fid, 'A:');
		for(cnt=1:length(Pos))
			fprintf(fid, '\t%f', A(cnt));
		end
		fprintf(fid, '\n');
		fprintf(fid, 'C:');
		for(cnt=1:length(Pos))
			fprintf(fid, '\t%f', C(cnt));
		end
		fprintf(fid, '\n');
		fprintf(fid, 'G:');
		for(cnt=1:length(Pos))
			fprintf(fid, '\t%f', G(cnt));
		end
		fprintf(fid, '\n');
		fprintf(fid, 'T:');
		for(cnt=1:length(Pos))
			fprintf(fid, '\t%f', T(cnt));
		end
		fprintf(fid, '\n');
		fclose(fid);
	end
end
