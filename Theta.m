function [theta] = Theta(L,sigma)
%
% calculates the spatial continuity penalty matrix theta according to the spatial coordinates of centroids of the anatomical regions for paper Mixtures of Block Models for Brain Networks 
% (By Zilong Bai, KDD Lab @ University of California, Davis)
%
% Input
%	L: N x d matrix. N row vectors, each being the spatial coordinate of the centroid of an anatomical region.
%	sigma: bandwidth for the penalty kernel calculation. It impacts the size of the blocks that each consists of spatially closely located anatomical regions.
% Output
%	theta: N x N penalty matrix for each pair of anatomical regions
% Please refer to Unsupervised Network Discovery for Brain Imaging Data, Z. Bai, et al. KDD 2017 for more about the significance of Spatial Continuity Regularization in learning interpretable blocks from brain networks without additional supervision.

X = L;
Y = L;

[xx,xy] = size(X);
[yx,yy] = size(Y);

XY = [];
for i = 1:xx
for j = 1:yx
  XY(i,j) = (norm(X(i,:)-Y(j,:)))^2;
end
end

theta = exp(XY./(2*sigma^2));

