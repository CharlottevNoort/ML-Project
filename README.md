# ML-Project

Scripts:

MLProjectCode_ALL.m
Main script for preprocessing and analysis of full dataset.
Takes as input full original matrix (data_breastcancer.mat).

MLProjectCode_PAM.m
Main script for preprocessing and analysis of PAM50 proteins.
Takes as input the recuced dataset that includes only these proteins (data_PAM50_77patients.mat).

mncn.m
Mean centering function.
Taken from Biosystems Data Analysis course (University of Amsterdam).

pcaKFold.m
Principal Component Analysis cross-validation code.

Data:

data_breastcancer.mat
Matrix containing original data, as extracted from .csv file obtained from kaggle.com.

data_PAM50_77patients.mat
Matrix conatining original data for the PAM50 proteins, as extracted from data_breastcancer.mat based on list of PAM50 genes and corresponding proteins.

PAM50_groups.mat
List of groups to which each patient belongs according to PAM50-mRNA analysis.
Extracted from additional .csv file on clinical data of the patients, obtained from kaggle.com.
