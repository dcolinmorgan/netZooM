% Set Program Parameters
% expression matrix and descriptors
exp_file   = 'gene_exp_v14.txt';
gene_file  = 'genes2.txt';
exp_subj   = 'EXP_patients.txt';
% generic motif file with CGs matched to TF-gene link
motif_file = 'MotifPriors/Motif3.txt';
% methylation array beta values with descriptors
methy_file = 'betas.clean_1436199620_V13.txt';
methyl_CGs = 'CGs.txt';
methyl_subj= 'CGpatients.txt';
% ppi and I/O paths
ppi_file   = 'ppi.txt';
panda_out  = '~/';  % optional, leave empty if file output is not required
save_temp  = '';  % optional, leave empty if temp data files are not needed afterward
lib_path   = '~/netZooM';  % path to the folder of PANDA source code
alpha      = 0.1;
