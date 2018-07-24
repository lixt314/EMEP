# EMEP

Single-Cell RNA-seq Data using Evolutionary Multiobjective Ensemble Pruning

In this study, a novel evolutionary multiobjective ensemble algorithm that imposes a specific structure on the ensemble clustering pruning, motivation by the observation that not all basic partition clustering results are very suitable for the ensemble clusters. In the algorithm, the unsupervised dimensionality reduction method is employed to project data from original high dimensional spaces to lower dimensional subspaces. After that, three different cluster validity indices including the overall cluster deviation, compactness, and the number of chosen basic partition clusters are proposed as objective functions to model multiple characteristics of the evolving clusters. After that, evolutionary multiobjective ensemble pruning algorithm is proposed to remove some clusters from the cluster set of an ensemble and further improve the generalization performance.

MainDE.m is the main function.

Objectivefucntion.m includes three objective function used in this work.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The six small-scale data sets are provided in the directory Data (https://github.com/ishspsy/project/tree/master/MPSSC).

Data_Deng.mat refers to http://science.sciencemag.org/content/343/6167/193.

Data_Ting.mat refers to https://www.ncbi.nlm.nih.gov/pubmed/25242334.

Data_Treutlin.mat refers to https://www.ncbi.nlm.nih.gov/pubmed/24739965.

Data_Ginhoux.mat refers to https://www.ncbi.nlm.nih.gov/pubmed/26054720.

Data_Buettner.mat refers to https://www.ncbi.nlm.nih.gov/pubmed/25599176.

Data_Pollen.mat refers to https://www.nature.com/articles/nbt.2967.

For the large scale data, Data_Zeisel.mat refers to https://github.com/BatzoglouLabSU/SIMLR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DOWNLOAD

We provide MATLAB implementations of EMEP.

Authors

Xiangtao Li and Ka-Chun Wong

Department of Information Science and Technology, Northeast Normal University, Changchun, Jilin, China
Department of Computer Science, City University of Hong Kong, Hong Kong

Contact

lixt314@nenu.edu.cn

