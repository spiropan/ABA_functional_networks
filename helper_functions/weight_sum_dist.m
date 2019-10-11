function [weighted_sums] = weight_sum_dis(W_edges_distances);
% function [X,G] = weight_sum_dis(W_edges_distance);
% Applies a weighted average to a vector of distances to create a model
% that will predict frequency of significant SFs

dis_bins=[0:8:160]; 
% Theses are the mean tissue-correlations for distances (shown in Figure
% 1). 
weights_vec=[0.1113,0.0983,0.0901,0.0853,0.0827,0.0783,0.0745,0.0720,0.0711,0.0688,...
            0.0679,0.0671,0.0664,0.0646,0.0653,0.0644,0.0653,0.0667,0.0630,0.0637];
        
for r=1:length(W_edges_distances)
    tmp=[];
    for d=1:length(dis_bins)-1
        inds=find(W_edges_distances{r}>dis_bins(d) & W_edges_distances{r}<dis_bins(d+1));
        tmp(d)=sum(W_edges_distances{r}(inds))*weights_vec(d);
    end
    weighted_sums(r)=sum(tmp);
end

