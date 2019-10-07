function [median_real, median_rand] = plot_distances(ind,ind_rand,D);
% Use this function to plot average distances of the real vs null networks
% as in Figure 1 of the 2017 reply by Richiardi

[X_real,G_real]=return_XG(ind,D);

% Create vector for the distances from the null networks
X_rand=[]; G_rand=[];
for n=1:size(ind_rand)
    [X_tmp,G_tmp]=return_XG(ind_rand(n),D);
    X_rand=[X_rand; X_tmp];
    G_rand=[G_rand; G_tmp];
end

X=[X_real; X_rand];
G1 = [G_real; G_rand]; 
G2 = [ones(length(G_real),1)*1; ones(length(G_rand),1)*2];
G=[G1 G2];
boxplot(X,G);

disp('median distance of the real Wi')
median_real=median(X_real)

disp('median distance of the null Wi')
median_rand=median(X_rand)