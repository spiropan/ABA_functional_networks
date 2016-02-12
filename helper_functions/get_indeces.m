function [ind]=get_indeces(numregions)
% Use this to quickly get lower triangle inds of symmetric matrix for N
% regions

oner=ones(numregions,numregions); oner=tril(oner,-1); ind=find(oner==1);