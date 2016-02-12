function []=spider_SVM_plot(input,x,labels,measure)
% function []=spider_SVM_plot(input,x,labels,measure)
% Use this to create plot of SVM accuracies plus null.
% Put multiple experiments (.mat file) into cell array to
% load them and plot together
% x is topn, labels is for legend and should follow order of cell array
% This was updated to now include 90% CI computed using the null values,
% computed non-parametrically
 
shapes={'.','s','*','d','o','v'}
lines={'-k','-k'}
specificity=0;

kill_figures;
meaner=[];
fontsize=16;
markersize=10;
alpha=0.1 % This to compute 90% CI

for e=1:length(input) % go over all elements of the struct
    load(input{e})
    if isfield(results,'sensitivity')
        results.bACC=(results.sensitivity+results.specificity)/2
    end
    % here also plot the chance line
    if isfield(results,'num_each_group')
        ylabeler='Classification';
        if specificity==1
            y=1-(1/length(results.num_each_group))*ones(1,length(x));
        else
            y=(1/length(results.num_each_group))*ones(1,length(x));
        end
    elseif isfield(results,'R')
        ylabeler='Prediction';
        y=[];
    else
        ylabeler='Classification';
        y=0.5*ones(1,length(x));
    end
    
    if nargin > 3
        results.accuracy=getfield(results,measure);
    else
        if isfield(results,'AUC')
            results.accuracy=results.AUC;
            measure='AUC'
        else
            measure='Accuracy'
        end
    end
    
    %[h,p,ci,stats] = ttest(input(e).AccNull,0.5,0.05);
    x
    if size(results.accuracy,1) > 1
        meaner=[meaner; results.accuracy(2:end,x)];
        plot_null=1;
    else 
        plot_null=0;
    end
    hold on; p=plot(x,results.accuracy(1,x),['-k' shapes{e}],'MarkerSize',markersize);
    %hold on; p=plot(x,results.accuracy(1,x),[lines{e}],'LineWidth',2);
end
if plot_null==1
    mean_null=mean(meaner);
    iter=size(meaner,1);
    % here compute nonParametric CI
    for f=1:size(mean_null,2)
        temp=sort(meaner(:,f),1,'ascend');
        low_val=temp(ceil(alpha/2*iter));
        hi_val=temp(ceil((1-alpha/2)*iter));
        ci_low(1,f)=abs(mean_null(1,f)-low_val);
        ci_hi(1,f)=abs(mean_null(1,f)-hi_val);
    end
    %std_null=std(meaner);
    %meaner=[]
    %hold on; errorbar(x,mean_null,ci(1,:),ci(2,:));
    %hold on; errorbar(x,mean_null,mean_null-std_null,mean_null+std_null);
    if isempty(meaner) ~= 1
        %hold on; errorbar(x,mean_null,std_null,std_null,'-k');
        hold on; errorbar(x,mean_null,ci_low,ci_hi,'-k');
        labels{end+1}='Null';
    end
end
%set(gca,'XTick',[0:5:x(end)])
%hold on; plot(x,y,'--k') % plot the 50% chance line
%set(p,'LineWidth',2)
%set(gca,'FontName','FixedWidth')

xlim([0 x(end)+1])
xlabel('Number of features','fontsize',fontsize)
ylabel([ylabeler ': ' measure],'fontsize',fontsize)
%title('Classifications','fontsize',fontsize);
legend(labels,'fontsize',fontsize,'Location','NorthEast')
