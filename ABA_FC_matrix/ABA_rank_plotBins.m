function [] = ABA_rank_plotBins(models,varargin)
% function [] = ABA_rank_plotBins(models,varargin)
% Here models is cell array of file name of .mat containing the ranked gene lsts
% Folder is the base folder containing the models (i.e. /nfs/zorba/INDI/TC)
% This will assume putting in 1 models, and plotting for all distance_bins 
% name is a string
% By defaut will plot for positive but change line below to load the neg
% instead
% pass 'file' which is .mat file containing the ranked gene list
% pass 'bdir' for the base directory (default /nfs/zorba/INDI/TC/)
% pass 'figname' for figure tiel
i=1;

% set defaults
bdir='/nfs/zorba/INDI/TC/'
figname='Mantel R vs. removed gene'

while i<=numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'bdir'}
                bdir = varargin{i+1};
                i=i+2;
            case {'figname'}
                figname=varargin{i+1};
                i=i+2;
        end
        
    end
end

for m=1:length(models)
    load([bdir models{m}])
    %load([bdir '/' models{m} '_RS_bins_ranks_pos'])
    %figname=models{m};
    %figname(end-3:3)=[];
    %figname=str_replace(figname,'_',' ');
    
    for d=1:size(Rcum,2)
        peak_ind=find(Rcum(:,d)==max(Rcum(:,d)));
        num_genes(d)=length(Rcum(:,d))-peak_ind;
    
        % Now plot the histogram
        if size(Rcum,2) > 1
            subplot(3,3,d)
        else
            subplot(1,1,m)
        end
        
        plot(Rcum(:,d));
     
        %ntitle({[' disBin: ' int2str(d)]; [' Max R=' num2str(max(Rcum(:,d)))]; [ ' at ' int2str(num_genes(d)) ' genes'];},'location','northwest')
        ntitle({[' disBin: 32-144 mm']; [' Max R=' num2str(max(Rcum(:,d)))]; [ ' at ' int2str(num_genes(d)) ' genes'];},'location','northwest')
        
        ylabel('Mantel R')
        xlabel('Removed gene index (ranked least to most informative)')
        %ylim([0,0.32])
        
    end
end

figtitle([figname])

