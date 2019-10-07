function [out] = random_clusters();
% This assumes Cortical_MNI_coords.csv file is in same directory as the
% function

% Set these parameters to roughly match the parameters of the real dataset
% (i.e. about 500 samples within networks, rest out of network, with about 12
% within nextwork clusters

% UPDATE: Richiardi et. al. in their reply noted that random clusters have
% shorter distances on average than resting fMRI networks. So we'll
% tweak some parameters to make the random networks more similar to resting
% fMRI networks. We'll do this by tripling the initial number of random
% clusters, then combine them to form 13 'networks'. Based on 
% Figure 1 in Richiardi et. al. 2015, each Wi network is comprised of 2-4 
% spatially distributed cluster. Therefore we'll combine 3 clusters to generate 
% each random network, and then compare and plot the average distances to
% resting state fMRI networks as in Figure 1 of the 2017 reply by
% Richiardi. 

% Set below parameters. Here each network is comprised of 3 clusters. 
% UPDATE: 
num_clusters=39;
radius=14;
% UPDATE:
num_networks=13;


dis_btw_clus=radius*2; % at least twice the radius to make sure points are non-overlapping
coords=load('Cortical_MNI_coords.csv');
num_pts=size(coords,1);

for c=1:num_clusters,
    % here check to make sure the new point is far enough 
    % from all previous points
    if c > 1
        disp('assigning new cluster')
        disp(int2str(c))
        dd=0; % initialize the distance between clus_centers
        while min(dd) <= dis_btw_clus % Keep finding new ind until distance is greater than dis_btw_clus
            ind=ceil(rand(1)*num_pts); % Find a random cluster center
            for p=1:size(clus_centers,1)                       
                dd(p)=euc_dis(coords(ind,:),clus_centers(p,:)); % here dd is vector of distances to previous clus_centers             
            end
        end
        clus_centers(c,:)=coords(ind,:);
    else
        disp('assigned first cluster')
        % This for the first point
        ind=ceil(rand(1)*num_pts);
        clus_centers(c,:)=coords(ind,:); 
    end
    
    % now find all points within a radius of this point
    cn=1;
    for x=1:size(coords,1)
        dd=euc_dis(clus_centers(c,:),coords(x,:));
        if dd<radius
            cluster{c}(cn)=x;
            cn=cn+1;
        end       
    end
end

% UPDATE
% Here combine the clusters to make larger 'networks' with longer distances
% between the comprising clusters
tmp=randperm(num_clusters);
starti=1;
for i=1:num_networks
    inds=tmp(starti:starti+2);
    super_cluster{i}=cat(2,cluster{[inds]});
    starti=starti+3;
end
cluster=super_cluster; % replace cluster with super_clusters

% Now get the rest of the indeces that aren't in the above
all_wi=[];
for c=1:length(cluster)
    % This to find the largest clusters
    all_wi=[all_wi cluster{c}];
end
[ZrestofBrain]=setdiff([1:size(coords,1)],all_wi);

% Here make the final cell array which will replace the 
% ind.W_all cell array in the real analysis. 
% ZrestofBrain should be the 9th array in the cell
out.W_all=[cluster(1:8) ZrestofBrain cluster(9:end)];

% Here define the within network clusters to be the largest four
for c=1:length(out.W_all)
    size_clus(c)=length(out.W_all{c});
end
[Y,I]=sort(size_clus,'descend');
out.Wi=I(2:5); % Start from two because ZrestofBrain is largest
out.W=[1:8,10:13]; % Keep the same as the real data

