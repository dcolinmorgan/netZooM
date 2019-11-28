function[]=milipeed_run(start,stop)
	%%screen qrsh -l matlab=1 -l h_vmem=100111M 
	% exp=textread('gene_exp_v14.txt');
	% clear all;clc;
	% clearvars -except beta
	addpath(genpath('~/netZooM'));alpha = 0.1;

	disp('Reading in ppi data!');
	[TF1, TF2, weight] = textread('ppi.txt', '%s%s%f');
	TFCoop = eye(length(unique(TF1)));
	TFNames=unique(TF1);
	[~,i] = ismember(TF1, TFNames);
	[~,j] = ismember(TF2, TFNames);
	TFCoop(sub2ind([size(TFNames,1), size(TFNames,1)], i, j)) = weight;
	TFCoop(sub2ind([size(TFNames,1), size(TFNames,1)], j, i)) = weight;
	TFCoop = NormalizeNetwork(TFCoop);
	fprintf('%d PPIs!\n', length(weight)); %%84388 PPIs!

	%% prepare motif info
	[TF, gene, CG]=textread('MotifPriors/Motif3.txt', '%s%s%s'); %%690 x 21528
	beta=textread('betas.clean_1436199620_V13.txt');
	% beta=zeros(349826,160);
	CGname=textread('CGs.txt', '%s');
	[f,loc]=ismember(CG, CGname);
	TF=TF(f); gene=gene(f); %% 690 x 20173 %%restrict RegNet to beta CG list

	beta=beta((loc(loc>0)),:); %% 569437

	%%Sort subjects
	exp_subj= textread('EXP_patients.txt', '%s');
	cg_subj= textread('CGpatients.txt', '%s');
	[ff,locc]=ismember(cg_subj,exp_subj); %%restrict beta subjects to exp list

	%% Build Coexpression matrix
	disp('Reading in expression data!'); 
	% [TF, gene, CG]=textread('MotifPriors/MotifPriorRefSeq_r10_p4_-750_250.txt', '%s%s%s');
	Genes= textread('genes2.txt', '%s');

	% disp(['NumGene= ',num2str(length(unique(gene)))]);

	[fff,loccc]=ismember((TF), unique(TF1)); %% from 690 to 644 TF, from 20173 to 14582 genes
	TF=TF(fff);gene=gene(fff);
	% disp(['NumGene= ',num2str(length(unique(gene)))]);

	[f,loc]=ismember(gene, unique(Genes)); %%restrict GeneCoReg to RegNet gene list from 20173 to 14587 genes
	TF=TF(f);gene=gene(f);
	% disp(['NumGene= ',num2str(length(unique(gene)))]);

	[w]=ismember(Genes, unique(gene)); %%restrict GeneCoReg to RegNet gene list

	Exp=textread('gene_exp_v14.txt'); %%23627 x 151
	Exp=Exp(w,:);
	[NumGenes, NumConditions] = size(Exp);
	fprintf('%d genes and %d conditions!\n', NumGenes, NumConditions); 
	Exp = Exp';
	% exp=exp';
	corex=corr(Exp); %15116

	for(a=start:stop)
	MILIPEED(Exp,beta,corex,a,loc,TF,gene,locc,fff,f,TFCoop,alpha,exp_subj);
	end

end


