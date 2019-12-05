function milipeed_run(START,END)
	%%screen qrsh -l matlab=1 -l h_vmem=100111M 
	% exp=textread('gene_exp_v14.txt');
	% clear all;clc;
    % rnd(12); %seed
	% clearvars -except beta
	addpath(genpath('~/netZooM'));alpha = 0.1;
	addpath(genpath('~/analyses/MILIPEED'));

%% ============================================================================
%% Set Program Parameters and Path
%% ============================================================================
% Run configuration to set parameter first (e.g., run('panda_config.m');)
fprintf('Input expression file: %s\n', exp_file);
fprintf('Input motif file: %s\n', motif_file);
fprintf('Input PPI file: %s\n', ppi_file);
fprintf('Output PANDA network: %s\n', panda_out);
fprintf('Output temp data: %s\n', save_temp);
fprintf('Alpha: %.2f\n', alpha);
addpath(lib_path);
%% ============================================================================
%% Read in Data
%% ============================================================================
disp('Reading in PPI data!');
tic

	[TF1, TF2, weight] = textread(ppi_file, '%s%s%f');
	TFCoop = eye(length(unique(TF1)));
	TFNames=unique(TF1);
	[~,i] = ismember(TF1, TFNames);
	[~,j] = ismember(TF2, TFNames);
	TFCoop(sub2ind([size(TFNames,1), size(TFNames,1)], i, j)) = weight;
	TFCoop(sub2ind([size(TFNames,1), size(TFNames,1)], j, i)) = weight;
	TFCoop = NormalizeNetwork(TFCoop);
	fprintf('%d PPIs!\n', length(weight)); %%84388 PPIs!
toc

disp('Reading in methylation data!');
tic
	%% prepare motif info
	[TF, gene, CG]=textread(motif_file, '%s%s%s'); %%690 x 21528
	beta=textread(methy_file);
	% beta=randn(349826,160);
	CGname=textread(methyl_CGs, '%s');
	[f,loc]=ismember(CG, CGname);
	TF=TF(f); gene=gene(f); %% 690 x 20173 %%restrict RegNet to beta CG list

	beta=beta((loc(loc>0)),:); %% 569437
toc 
%%Sort subjects
exp_subj= textread(exp_subj, '%s');
cg_subj= textread(methyl_subj, '%s');
[ff,locc]=ismember(cg_subj,exp_subj); %%restrict beta subjects to exp list
locc(locc==0) = [];

%% Build Coexpression matrix
disp('Reading in expression data & building aggregate coexpression network!'); 
tic
	% [TF, gene, CG]=textread('MotifPriors/MotifPriorRefSeq_r10_p4_-750_250.txt', '%s%s%s');
	Genes= textread(gene_file, '%s');

	% disp(['NumGene= ',num2str(length(unique(gene)))]);

	[fff,loccc]=ismember((TF), unique(TF1)); %% from 690 to 644 TF, from 20173 to 14582 genes
	TF=TF(fff);gene=gene(fff);
	% disp(['NumGene= ',num2str(length(unique(gene)))]);

	[f,loc]=ismember(gene, unique(Genes)); %%restrict GeneCoReg to RegNet gene list from 20173 to 14587 genes
	TF=TF(f);gene=gene(f);
	% disp(['NumGene= ',num2str(length(unique(gene)))]);

	[w]=ismember(Genes, unique(gene)); %%restrict GeneCoReg to RegNet gene list

	Exp=textread(exp_file); %%23627 x 151
    % Exp=randn(23627,151);
	Exp=Exp(w,:);
	[NumGenes, NumConditions] = size(Exp);
	fprintf('%d genes and %d conditions!\n', NumGenes, NumConditions); 
	Exp = Exp';
	% exp=exp';
	corex=corr(Exp); %15116
toc
%% ============================================================================
%% Run MILIPEED
%% ============================================================================

	for(a=START:END)
	MILIPEED(Exp,beta,corex,a,loc,TF,gene,locc,fff,f,TFCoop,alpha,exp_subj);
	end

end


