function [out] = plot_correlations_cluster_size(files,labels,files_dis)
% function [out] = plot_correlations_cluster_size(files)
% files          = cell array of .mat files that are saved in the last cell
%                  of proc_all.m
% labels         = cell array to label the plots in each file (i.e. 'whole brain','omit
% Z rest of brain' etc.
% Pass files_dis which are runs with few null_size by 100 random networks
% (useful for quickly generating edge distance metrics

if nargin < 3
    files_dis=[];
end

markers={'ko-','b+-','r*-'}
line_markers={'k-','b-','r-'}
real_Wi=201;
real_W=501;
real_median_dist=50.97;
kill_figures

labels=[labels 'rsfMRI'];

h=gobjects(length(files)+1,1);

for f=1:length(files)
    
    if ~isempty(files_dis),
        tmp=load(files_dis{f});
    else
        tmp=load(files{f});
    end
    
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
        title('Simulated network properties vs. contiguous cluster size')
        legend(labels,'location','northwest')
    end
    
    % ------------ MIDDLE PANEL ----------------------------------
    subplot(3,1,2); hold on;
    h(f)=plot(tmp.clus_radius,tmp.median_dist_rand,markers{f});
    [means,ci_low,ci_hi]=calc_ci(tmp.medians_dist_rand);
    
    
    if f==length(files)
        % Here add the dashed line for the real rsfMRI network
        hold on; x=[5:0.1:16];
        h(f+1)=plot(x,ones(length(x),1)*real_median_dist,'--k','LineWidth',0.5)
        
        ylabel('Median distance (mm)')
        xlim([5,16])
        ylim([0,150])
        title('')
        hold on; e=errorbar(tmp.clus_radius,means,ci_low,ci_hi,'-k');
        legend(h,labels,'location','northwest')
        %legend(labels,'location','northwest')
    end
    
    % ------------ BOTTOM PANEL ----------------------------------
    % here can't use the dis (one shot) files
    tmp=load(files{f});
    
    subplot(3,1,3); hold on;
    plot(tmp.clus_radius,tmp.perc_sig*100,markers{f});
    
    if f==length(files)
        xlabel('Contiguous cluster radius (mm)')
        ylabel('% significant SFs')
        xlim([5,16])
        ylim([0,100])
        title('')
        
        legend(labels(1:end-1),'location','northwest')
    end
end

end

function [means,ci_low,ci_hi] = calc_ci(input_cell)
    alpha=0.1 % This to compute 90% CI

    for f=1:length(input_cell)
        iter=length(input_cell{f});
        temp=sort(input_cell{f},1,'ascend');
        low_val=temp(ceil(alpha/2*iter));
        hi_val=temp(ceil((1-alpha/2)*iter));
        ci_low(1,f)=abs(mean(input_cell{f})-low_val);
        ci_hi(1,f)=abs(mean(input_cell{f})-hi_val);
        means(1,f)=mean(input_cell{f});
    end
end
