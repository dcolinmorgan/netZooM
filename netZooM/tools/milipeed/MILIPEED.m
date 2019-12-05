function[]=MILIPEED(Exp,beta,corex,iii,loc,TF,gene,locc,fff,f,TFCoop,alpha,exp_subj)
%1:size(Exp,1))	
% iii=1;
	disp(['Computing coexpression network for subject ',num2str(iii), ':']);
	GeneCoReg=corr(Exp([1:iii-1 iii+1:size(Exp,1)],:));
	GeneCoReg= 1/size(Exp,1)*(corex- GeneCoReg )+GeneCoReg;
	GeneCoReg = NormalizeNetwork(GeneCoReg);				%%%14587 x 14587
	% dlmwrite(['corLP_',char(exp_subj(iii)), '.txt'],z)
	beta_pat=beta(:,locc(iii));
	Wcase=1-beta_pat;
	Wcase=Wcase(fff);Wcase=Wcase(f);
	% Wcase=beta_pat((loc(loc>0))); %% 569437
	% Wcase=Wcase(f); %% 543052


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
	GeneNames=unique(genenet);
	NumTFs = length(TFNames);
	[~,i] = ismember(TFnet, TFNames);
	[~,j] = ismember(genenet, uGene);
	RegNet = zeros(NumTFs, length(uGene));
	RegNet(sub2ind([NumTFs, length(uGene)], i, j)) = weight;
	fprintf('%d TFs and %d edges!\n', NumTFs, length(weight));
	RegNet = NormalizeNetwork(RegNet);
	disp('Running PANDA algorithm:');
	AgNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha);
	% dlmwrite([char(exp_subj(iii)),'_milipeed.txt'],AgNet(:))
	% SavePairs(TFNames, uGene, AgNet, RegNet, char(exp_subj(iii)))
	% SavePairs(TFNames, uGene, AgNet, RegNet, char(strcat('zero_',exp_subj(iii))));
	SavePairs(TFNames, uGene, AgNet, RegNet, char(fullfile('eye_MILIPEED',strcat('eye_',exp_subj(iii)))));

end