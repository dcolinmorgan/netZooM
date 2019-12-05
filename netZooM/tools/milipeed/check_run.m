exp_subj= textread('EXP_patients.txt', '%s');

a=dir('*.pairs')
for i=1:length(a)
	b=strsplit(a(i).name,'_');
	c(i)=b(2);
end
[d,e]=setdiff(exp_subj,c)