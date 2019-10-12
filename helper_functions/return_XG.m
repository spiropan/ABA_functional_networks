function [X,G] = return_XG(ind,D);
% function [X,G] = return_XG(ind,D);
% Creates the X and G (grouping) variable for inputs to boxplot (Figure 2)
% initialize X and G, then iterate over each of 13 networks to plot
G=[]; X=[];
nets_ind=ind.W;
for n=1:length(nets_ind)
    N=ind.W_all{nets_ind(n)}; % store the indeces of the nodes for each network
    net_dis=D(N,N); % isolate the distances for nodes within each network
    % extract lower triangle of those distances and turn into vec
    tmp=get_indeces(length(N));
    dis_vec=net_dis(tmp);

    % Create the X and G (grouping) variable for inputs to the boxplot
    X=[X; dis_vec];
    G=[G; ones(length(dis_vec),1)*nets_ind(n)];
end
