[TF, gene, CG]=textread('MotifPriors/MotifPriorRefSeq_mQTL_p4_-750_250.txt', '%s%s%s');
% tag='1';
% [CGname,Bcase,Bcont]=textread('../SPIDER/MedianBetaValues.txt', '%s%f%f', 'delimiter', '\t', 'headerlines', 1);
% exp=textread('exp.clean_V14_1435701821.txt');
% genes= textread('genes.txt', '%s');
% exp_subj= textread('EXP_patients.txt', '%s');

beta=textread('betas.clean_1436199620_V13.txt');
CGname=textread('CGs.txt', '%s');
CGpat=textread('CGpatients.txt', '%s');

[f,loc]=ismember(CG, CGname);
TF=TF(f); gene=gene(f);

for(ii=1:size(beta,2))
	beta_pat=beta(:,ii);
	beta_pat=1-beta_pat;

	Wcase=beta_pat(loc(loc>0));
% Wcont=Bcont(loc(loc>0));C
	
	% z=1/151(np.corrcoef(exp)- np.corrcoef(exp[:,0:150]))+np.corrcoef(exp[:,0:150])
	
	uGene=unique(gene);
	uTF=unique(TF);
	[~,i]=ismember(TF,uTF);
	[~,j]=ismember(gene, uGene);
	[net,~,idx]=unique([i,j], 'rows');
	NumEdges=size(net,1);
	muWcase = accumarray(idx,Wcase,[],@mean);
% muWcont = accumarray(idx,Wcont,[],@mean);

fid=fopen(['mili_',char(CGpat(ii)), '.pairs'], 'wt');
	for(ecnt=1:NumEdges)
		fprintf(fid, '%s\t%s\t%f\n', uTF{net(ecnt,1)}, uGene{net(ecnt,2)}, muWcase(ecnt));
	end
	fclose(fid);
end

% fid=fopen(['SPIDERmeth', tag, '_zero_controls.pairs'], 'wt');
% for(ecnt=1:NumEdges)
% 	fprintf(fid, '%s\t%s\t%f\n', uTF{net(ecnt,1)}, uGene{net(ecnt,2)}, muWcont(ecnt));
% end
% fclose(fid);
