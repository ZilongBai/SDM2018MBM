# SDM2018MBM Mixtures of Block Models for Brain Networks
This repository contains matlab source code to solve the Block Structure Set Linked Matrix tri-Factorization with Spatial Continuity Regularization (BSSLMtF-SCR) to learn sets of weighted mixing matrices M and small set of multiple block structures F (the main function for implementing the iternative solver framework is mbmsolver_continuity.m). Each block structure can be regularized with Spacial Continuity depending on the application domain (e.g. Brain Imaging Data). Through rescaling function (i.e., AlphaExtractor_normM.m), the weights for each block model in modeling each single graph that comprise the multigraph X can be extracted, which achieves graph embedding.

The details of the model and the solver are described in our paper [Mixtures of Block Models for Brain Imaging Data](https://epubs.siam.org/doi/10.1137/1.9781611975321.6), which demonstrates their utility in the analyses of both synthetic networks constructed with randomly generated underlying block structures, mixing matrices, and block models weights, as well as networks constructed based on real-world brain scans. 

The details of deductions for our method as well as the techniques from related work that are involved in developing our method can also be found in our paper. See the references of our paper for more information on related work.

Please refer to [Unsupervised Network Discovery for Brain Imaging Data](http://dl.acm.org/citation.cfm?id=3098023&CFID=796408940&CFTOKEN=92880021) for more about the significance of Spatial Continuity Regularization in learning interpretable blocks from brain networks without additional supervision.

File(s) in this repository:

updateXbar.m : serves to calculate intermediate variables, in particular graphs reconstructed with current F and M, to reduce computationl redundancy in implementing update rules for F and M. 

Theta.m : calculates the spatial continuity penalty matrix theta according to the spatial coordinates of centroids of the anatomical regions for paper Mixtures of Block Models for Brain Networks.

AlphaExtractor_normM.m : extracts alphas from weighted mixing matrices M (Equation 4.17 of the paper Mixtures of Block Models for Brain Networks, SDM 2018).       

reconserrs.m : serves to calculate the relative reconstruction error of the multigraph.      

updateF_continuity.m : serves to use multiplicative update rules to update nonnegative tensor F (Equation 4.15 in paper Mixtures of Block Models for Brain Networks, SDM 2018).

InitFMAlpha.m : randomly initializes F, M, and alpha given N, k, p, and r.     

updateM.m : serves to use multiplicative update rules to update nonnegative tensor M (equation 4.17 in paper Mixtures of Block Models for Brain Networks, SDM 2018).    

mbmsolver_continuity.m : solves the Block Structure Set Linked Matrix tri-Factorization with Spatial Continuity Regularization (BSSLMtF-SCR) for sets of  weighted mixing matrices M and small set of block structures F (for Algorithm 1 Mixtures of Block Models Discovery in Paper Mixtures of Block Models for Brain Networks, SDM 2018).

Please report any bugs/issues to zlbai@ucdavis.edu
