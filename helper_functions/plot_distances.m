% function [] = plot_avg_distances(real_ind,null_ind,D);
% Use this function to plot average distances of the real vs null networks
% as in Figure 1 of the 2017 reply by Richiardi

[X_real,G_real]=return_XG(ind,D);
[X_rand,G_rand]=return_XG(ind_rand,D);

X=[X_real; X_rand];
G1 = [G_real; G_rand] 
G2 = [ones(length(G_real),1)*1; ones(length(G_rand),1)*2]
G=[G1 G2];
boxplot(X,G);
