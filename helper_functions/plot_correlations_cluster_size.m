function [out] = plot_correlations_cluster_size(files,labels)
% function [out] = plot_correlations_cluster_size(files)
% files          = cell array of .mat files that are saved in the last cell
%                  of proc_all.m
% labels         = cell array to label the plots in each file (i.e. 'whole brain','omit
% Z rest of brain' etc.

markers={'ko-','g+-','b*-'}
line_markers={'k-','g-','b-'}
real_Wi=201;
real_W=501;
real_median_dist=50.97;
kill_figures

labels=[labels 'rsfMRI'];

for f=1:length(files)
    tmp=load(files{f});
    
    % ------------ TOP PANEL ------------------------------------
    subplot(3,1,1); hold on;
    %plot(tmp.clus_radius,tmp.median_W_samples,line_markers{f},tmp.clus_radius,tmp.median_W_samples,markers{f});
    plot(tmp.clus_radius,tmp.median_W_samples,markers{f});

    % Here add the dashed line for the real rsfMRI network
    if f==length(files)
        hold on; x=[5:0.1:16];
        plot(x,ones(length(x),1)*real_W,'--k','LineWidth',0.5)
        ylabel('Median # W nodes')
        xlim([5,16])
        title('')
        legend(labels,'location','northwest')
    end
    % ------------ MIDDLE PANEL ----------------------------------
    subplot(3,1,2); hold on;
    plot(tmp.clus_radius,tmp.median_dist_rand,markers{f});
    
    if f==length(files)
        % Here add the dashed line for the real rsfMRI network
        hold on; x=[5:0.1:16];
        plot(x,ones(length(x),1)*real_median_dist,'--k','LineWidth',0.5)

        ylabel('Median distance (mm)')
        xlim([5,16])
        ylim([0,150])
        title('')
        legend(labels,'location','northwest')
    end
    
    % ------------ BOTTOM PANEL ----------------------------------
    subplot(3,1,3); hold on;
    plot(tmp.clus_radius,tmp.perc_sig,markers{f});
    
    if f==length(files)
        xlabel('Simulated networks: cluster radius (mm)')
        ylabel('% significant SFs')
        xlim([5,16])
        ylim([0,1])
        title('')
        
        legend(labels(1:end-1),'location','northwest')
    end
end
  