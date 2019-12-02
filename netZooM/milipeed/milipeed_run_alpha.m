%%screen qrsh -l matlab=1 -l h_vmem=100111M 
% exp=textread('gene_exp_v14.txt');
clearvars -except beta
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
% beta=textread('betas.clean_1436199620_V13.txt');
CGname=textread('CGs.txt', '%s');
[f,loc]=ismember(CG, CGname);
TF=TF(f); gene=gene(f); %% 690 x 20173 %%restrict RegNet to beta CG list

[fff,loccc]=ismember((TF), unique(TF1)); %% 349 x 19796
TF=TF(fff); gene=gene(f);  %%restrict RegNet TF to beta ppi
%%Sort subjects
exp_subj= textread('EXP_patients.txt', '%s');
cg_subj= textread('CGpatients.txt', '%s');
[ff,locc]=ismember(cg_subj,exp_subj); %%restrict beta subjects to exp list

%% Build Coexpression matrix
disp('Reading in expression data!'); 
% [TF, gene, CG]=textread('MotifPriors/MotifPriorRefSeq_r10_p4_-750_250.txt', '%s%s%s');
Genes= textread('genes2.txt', '%s');
[f]=ismember(Genes, unique(gene)); %%restrict GeneCoReg to RegNet gene list
% TF=TF(f);gene=gene(f);

Exp=textread('gene_exp_v14.txt'); %%23627 x 151
Exp=Exp(f,:);
beta2=beta(fff,:);					%%% 14587 x 151
% Genes=
[NumGenes, NumConditions] = size(Exp);
fprintf('%d genes and %d conditions!\n', NumGenes, NumConditions); 
Exp = Exp';
% exp=exp';
corex=corr(Exp); %15116


[f]=ismember(gene, unique(Genes)); %%restrict GeneCoReg to RegNet gene list
TF=TF(f);gene=gene(f);

% function[]=MILIPEED(Exp,beta,corex,a,b,loc,TF,locc)
% for(iii=a:b)%1:size(Exp,1))	
iii=1;
	disp(['Computing coexpression network for subject ',num2str(iii), ':']);
	GeneCoReg=corr(Exp([1:iii-1 iii+1:size(Exp,1)],:));
	GeneCoReg= 1/size(Exp,1)*(corex- GeneCoReg )+GeneCoReg;
	GeneCoReg = NormalizeNetwork(GeneCoReg);				%%%14587 x 14587
	% dlmwrite(['corLP_',char(exp_subj(iii)), '.txt'],z)
	beta_pat=beta(:,locc(iii));
	beta_pat=1-beta_pat;

	Wcase=beta_pat((loc(loc>0))); %% 569437
	% Wcase=Wcase(f); %% 543052

	[fff,loccc]=ismember((TF), unique(TF1)); %% 349 x 19796
	TF=TF(fff);gene=gene(fff);Wcase=Wcase(fff);
	[f]=ismember(gene, unique(Genes)); %%restrict GeneCoReg to RegNet gene list
	TF=TF(f);gene=gene(f);Wcase=Wcase(f);

	uGene=unique(gene);  %% 569437 to 20173
	uTF=unique(TF);      %% 569437 to 690
	[~,i]=ismember(TF,uTF);
	[~,j]=ismember(gene, uGene);
	[net,~,idx]=unique([i,j], 'rows'); %% 442205x2
	NumEdges=size(net,1);
	disp(['Computing methylation motif network for subject ',num2str(iii), ':']);
	muWcase = accumarray(idx,Wcase,[],@mean);

	for(ecnt=1:NumEdges)
		TFnet{ecnt}=uTF{net(ecnt,1)};
		genenet{ecnt}=uGene{net(ecnt,2)};
		weight(ecnt)=muWcase(ecnt);
	end

	TFNames = unique(TFnet);
	NumTFs = length(TFNames);
	[~,i] = ismember(TFnet, TFNames);
	[~,j] = ismember(genenet, uGene);
	RegNet = zeros(NumTFs, length(uGene));
	RegNet(sub2ind([NumTFs, length(uGene)], i, j)) = weight;
	fprintf('%d TFs and %d edges!\n', NumTFs, length(weight));
	RegNet = NormalizeNetwork(RegNet);
	disp('Running PANDA algorithm:');
	AgNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha);
	dlmwrite([exp_subj(iii),'_milipeed.txt'],AgNet)

% end

