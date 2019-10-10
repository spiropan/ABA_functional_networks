%% Here set variables
% Below is directory where all ABA MA samples have been downloaded and
% unzipped
clear all
collect_data=0; % If already assembled MA data matrix keep this 0 to save time
adjust_brain_batch=0; % If already regressed out brain and batch id from MA matrix to save time
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';

addpath([files_dir 'helper_functions'])
% Import the date file for which samples were used in the science paper
% (S1). data(:,1) will have donor ids, dat(:,3) will have well id
[dat(:,1),dat(:,2),dat(:,3),dat(:,4),dat(:,5),dat(:,6),dat(:,7),dat(:,8),dat(:,9),...
    dat(:,10),dat(:,11),dat(:,12),dat(:,13)]=textread([files_dir 'Data_File_S1.txt'],'%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t');

% Import the date for for which genes were used in the Science paper (S2)
% gene, gen(:,1) will have the probe ids used to extract the right  expression
% values
[gen(:,1),gen(:,2),gen(:,3)]=textread([files_dir 'Data_File_S2.txt'],'%s%s%s','delimiter','\t');

if collect_data==1
    % Here extract the genes from each model
    models={'H0351_2001','H0351_2002','H0351_1009','H0351_1012','H0351_1015','H0351_1016'}; % in same order as in S1
    MA=[]; % this will initialize the microarray data matrix used in the study
    probes_tocheck=[]; samples_tocheck=[];
    
    for i=1:6
        disp(['working on donor: ' models{i}])
        model_wperiod=models{i}; model_wperiod(6)='.'; % Here replace '_' with a '.'
        
        % Here determine indeces of which samples to include from each donor
        % Extract the SampleAnnotation file
        sa=importdata([ABA_dir models{i} '/SampleAnnot.csv']);
        wellid=sa.textdata(2:end,3); % Extract all  wellids for this donor
        model_wperiod;
        subind=cell_exp_ind(dat(:,1),model_wperiod,0); % get indeces for each donor from master spreadsheet
        wellid_filt=dat(subind,2); % Here extract the wellids for this donor
        [c,ia,ib_samples]=intersect(wellid_filt,wellid,'stable');
        
        % Here extract the ~17,000 probes used in the paper, probes(:,2) has the
        % right ids
        [probes(:,1),probes(:,2),probes(:,3),probes(:,4)]=textread([ABA_dir models{i} '/Probes.csv'],'%s%s%s%s','delimiter',',');
        
        % Here determine indeces of probes to extract from MA data
        [c,ia,ib_probes]=intersect(gen(:,1),probes(:,2),'stable');
        
        tmp=importdata([ABA_dir models{i} '/MicroarrayExpression.csv']);
        tmp(:,1)=[]; % Get rid of the 1st column which are numerical probe ids
        
        % Here the row indeces select out the right probes/genes, and the column
        % indeces will select out the 1777 samples used in the Richiardi study (when
        % combined across subjects)
        MA=[MA tmp(ib_probes,ib_samples)]; % here still have to determine the right ind variable to use
        probes_tocheck=[probes_tocheck; probes(ib_probes,2)]; % This to check with S2 file
        samples_tocheck=[samples_tocheck; wellid(ib_samples)];
    end
    save([ABA_dir '/Science_Paper_MA.mat'],'MA','dat','gen','ABA_dir','files_dir');
else
    load([ABA_dir '/Science_Paper_MA.mat']);
end

%% At this point MA is the total data matrix in same sample and probe order as in S1 and S2.
% Now regress out brain and and batchid as in the Science paper, brain id is
% in dat(2:end,1), batch_id in dat(2:end,3)
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA')
    load([ABA_dir '/Science_Paper_MA.mat']);
end
addpath([files_dir 'helper_functions'])

% Create regressor for brain id
adjust_brain_batch=1
if adjust_brain_batch==1
    brain_ids=unique(dat(2:end,1)); X1=zeros(size(dat(2:end,1),1),1); % Initialize brainid regressor
    for i=1:length(brain_ids)
        disp(['working on brain id: ' brain_ids{i}])
        donor_ind{i}=cell_exp_ind(dat(2:end,1),brain_ids{i},1);
        X1(donor_ind{i}',1)=i;
    end
    % Donor_ind{i} will come in handy for permutation procedure later
    
    % Create regressor for the batchid
    batch_ids=unique(dat(2:end,3)); X2=zeros(size(dat(2:end,3),1),1); %Initialize batchid regressor
    for i=1:length(batch_ids)
        disp(['working on batch id: ' batch_ids{i}])
        batch_ind{i}=cell_exp_ind(dat(2:end,3),batch_ids{i},1);
        X2(batch_ind{i}',1)=i;
    end
    X1D=dummyvar(X1); X2D=dummyvar(X2);
    
    % Now go over each gene (~17K) and regress out brain and batchid 
    for g=1:size(MA,1) % Go over all genes
        disp(['regressing brainid and batchid from gene: ' int2str(g)]) 
        % Here do brainid first and then batchid
        [b,dev,stats1]=glmfit(X1D(:,2:end),MA(g,:)');
        [b,dev,stats]=glmfit(X2D(:,2:end),stats1.resid);
        MA_resid(g,:)=stats.resid';
    end
    %save([ABA_dir '/Science_Paper_MA_resid.mat'],'MA_resid','donor_ind','dat','gen','ABA_dir','files_dir')
else
    load([ABA_dir '/Science_Paper_MA_resid.mat'])
end

%% Extract the indeces for 13 networks, 4 networks of interest and from each subject
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA_resid')
    load([ABA_dir '/Science_Paper_MA.mat']);
end
addpath([files_dir 'helper_functions'])

nets=unique(dat(2:end,13));
ind.donor_ind=donor_ind; % This to pass into below funciton for permuting within subjects

% In the above nets cell array, dDMN has index =10, Salience=6, Sensorimotor=7, Visuospatial=8,
% Zrest_of_brain=9
for n=1:length(nets)
    ind.W_all{n}=cell_exp_ind(dat(2:end,13),nets{n},1);
end
% Index the indeces of the above cell array
ind.Wi=[10,6,7,8];
ind.W=[1:8,10:length(nets)];

% Remove edges within same tissue class and create a "censor" matrix with
% ones to indicate tissue-tissue edges. Tissues classes are listed in dat(:,8)
tissue_classes=unique(dat(2:end,8)); censor_mat=zeros(size(MA,2),size(MA,2)); % This should be 1,777x1,777
for t=1:length(tissue_classes)
    disp(['removing edges for tissue: ' tissue_classes{t}])
    tind=cell_exp_ind(dat(2:end,8),tissue_classes{t},1);
    if length(tind) > 1
        pairs=nchoosek(tind,2); % This to generate all unique pairs among these indeces
        for p=1:size(pairs,1)
            censor_mat(pairs(p,1),pairs(p,2))=1; censor_mat(pairs(p,2),pairs(p,1))=1;
        end
    end
end

% Replicate primary analyses in Richiardi et. al. 
% Compute the total tissue similarity matrix, zero out negative edges and
% within-tissue edges
T_mat=corr(MA_resid); T_mat(find(T_mat<0))=0; 
T_mat(find(censor_mat==1))=0;
[results_resid]=compute_SF(T_mat,ind,10);
 %save([ABA_dir '/Science_Paper_MA_resid.mat'],'MA_resid','donor_ind','dat','gen','ABA_dir','files_dir',...
     %'censor_mat','results_resid')
 
%% Here compute distances between all samples to compute distance matrix.
% columns 10-12 have MNI coordinates, 13 has network name. Panel C in
% writeup
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA_resid')
    load([ABA_dir '/Science_Paper_MA_resid.mat']);
    load([ABA_dir '/Science_Paper_MA.mat']);
end
addpath([files_dir 'helper_functions'])

coords=[];
for i=1:size(dat(2:end,10)) % Go over all MNI x coordinates
    coords(i,1)=str2num(dat{i+1,10}); % X coord
    coords(i,2)=str2num(dat{i+1,11}); % y coord
    coords(i,3)=str2num(dat{i+1,12}); % z coord
end

D=zeros(length(coords),length(coords)); % Initialize adj matrix to hold distances
for i=1:size(D,1)
    for j=1:size(D,2)
        D(i,j)=sqrt((coords(i,1)-coords(j,1))^2+(coords(i,2)-coords(j,2))^2+(coords(i,3)-coords(j,3))^2);
        D(j,i)=D(i,j); % Just to make symmetric
    end
end

% Here remove edges according to min distances specificied in Dthr
Dthr=[4:4:24];
T_mat=corr(MA);
T_mat(find(T_mat<0))=0; % Here zero out negative edges
T_mat(find(censor_mat==1))=0; % Here remove the tissue-tissue-edges too

for d=1:length(Dthr)
    NewT_mat=T_mat;
    NewT_mat(find(D<=Dthr(d)))=0; % Here zero out distances <X mm.
    [results_MinDis(d).out]=compute_SF(NewT_mat,ind,200);
end

% Here plot the strength fractions for various corrections
results_plot(:,1)=[results_resid.real_SF; results_resid.null_SF']; % First add the Richiardi results
for d=1:length(Dthr)
    results_plot(:,d+1)=[results_MinDis(d).out.real_SF; results_MinDis(d).out.null_SF'];
end
results.accuracy=results_plot;
save('to_plot','results');
spider_SVM_plot({'to_plot.mat'},[1:7],{'Real'})
h_legend=legend('Real SF','Null SF')
set(h_legend,'fontsize',14)
h=gca;
set(h,'XTickLabel',{'','Tissue','<4 mm','<8 mm','<12 mm','<16 mm','<20 mm','<24 mm'},'fontsize',12)

%% Plot distance vs. tissue-tissue correlationin the within-network (Wi) and 
% out of network (T-W) edges. Panel B in writeup
% Wi          = dark grey, 
% T-W         = light grey, 
% Order in ind.Wi is dDMN, Salience, Sensorimotor, Visuospatial
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA_resid')
    load([ABA_dir '/Science_Paper_MA_resid.mat']);
end
addpath([files_dir 'helper_functions'])

T_mat=corr(MA_resid); N=size(T_mat,1);
D_mat=D; 

% Here zero out negative edges
D_mat(find(T_mat<0))=NaN; % Make sure to do D_mat first here
T_mat(find(T_mat<0))=NaN; 

% Here create a figure 
% here remove tissue-tissue connections
kill_figures
T_mat(find(censor_mat==1))=NaN;
D_mat(find(censor_mat==1))=NaN;

% here grab T_W indeces
disp('Working on T-W')
% First define T_vec (total vector of edges strengths and distances
T_ind=get_indeces(N);
T_vec=T_mat(T_ind); D_vec=D_mat(T_ind);

% Now collect vector indeces of all W within network edges to subtract
% from the T_vec
for i=1:length(ind.W)
    disp(['working on ' int2str(ind.W(i))])
    pairs=nchoosek(ind.W_all{ind.W(i)},2); % This to generate all unique pairs among these indeces
    pairs=[pairs; pairs(:,2) pairs(:,1)]; % here to add symmetric columns   
    % here convert the within network i,j edges to vector indeces
    [ind_to_remove]=sub2ind(size(T_mat),pairs(:,1),pairs(:,2));
    [c,ia,ib]=intersect(T_ind,ind_to_remove); % Here to find only the lower triangle vectorized indeces
    T_ind(ia)=[];
end

TW_vec=T_mat(T_ind); DTW_vec=D_mat(T_ind);
%scatter(DTW_vec, TW_vec,'bo');
scatter(DTW_vec, TW_vec,[],[0.7,0.7,0.7],'filled');

% Here grab all the Wi edges and distances
disp('working on Wi')
Wi_vec=[]; DWi_vec=[];
for i=1:length(ind.Wi)
    disp(['working on ' int2str(i)])
    % here extract the within network edges
    Wi{i}=T_mat(ind.W_all{ind.Wi(i)},ind.W_all{ind.Wi(i)});
    % here create vector of all 4 within network edge strengths 
    Wi_vec=[Wi_vec; Wi{i}(get_indeces(size(Wi{i},1)))];
    % Do the same for the distances    
    DWi{i}=D_mat(ind.W_all{ind.Wi(i)},ind.W_all{ind.Wi(i)});
    DWi_vec=[DWi_vec; DWi{i}(get_indeces(size(Wi{i},1)))];
end

%hold on; scatter(DWi_vec,Wi_vec,'go');
hold on; scatter(DWi_vec, Wi_vec,[],[0.5,0.5,0.5],'filled');

% Here add the best fit line
DWi_vec(find(isnan(DWi_vec)==1))=[];
Wi_vec(find(isnan(Wi_vec)==1))=[];
DTW_vec(find(isnan(DTW_vec)==1))=[];
TW_vec(find(isnan(TW_vec)==1))=[];

tot_exp_vec=[Wi_vec; TW_vec]; tot_dis_vec=[DWi_vec; DTW_vec];

% here plot best fit line
hold on; plot_linear(tot_dis_vec,tot_exp_vec,'k');
ylim([0,0.67])
xlim([0,180])
h_legend=legend('T-W','Wi');
set(h_legend,'fontsize',13);
xlabel('Euclidean Distance: mm','fontsize',16)
ylabel('Tissue-tissue correlation','fontsize',16)

% Here compute linear correlation and p-value
[B,dev,stats]=glmfit(tot_dis_vec,tot_exp_vec);
[rho,pval]=corr(tot_dis_vec,tot_exp_vec);

% here examine differences in distance between Wi and T-W
DWi_mean=mean(DWi_vec)
DTW_mean=mean(DTW_vec)
[h,p,ci,stats]=ttest2(DWi_vec, DTW_vec)

%% Here can try regressing out distance as well
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA_resid')
    load([ABA_dir '/Science_Paper_MA_resid.mat']);
end
addpath([files_dir 'helper_functions'])

T_mat=corr(MA_resid); N=size(T_mat,1);

% Here zero out negative edges and remove within-tissue edges
T_mat(find(T_mat<0))=NaN; T_mat(find(censor_mat==1))=NaN;
D_mat=D; 

% initialize matrix that will hold adjusted correlations
T_mat_dis=zeros(size(T_mat)); 

disp('Ok regressing out distance')
exp_vec=T_mat(get_indeces(N)); zero_out=find(isnan(exp_vec)==1); exp_vec(zero_out)=NaN;
dis_vec=D(get_indeces(N)); dis_vec(zero_out)=NaN;

[B,dev,STATS]=glmfit(dis_vec,exp_vec);
yresid=STATS.resid;
% Here use polyfit to fit cubic fit to the data
%p=polyfit(dis_vec,exp_vec,1);
%yfit=polyval(p,dis_vec);
%yresid=exp_vec-yfit;

% Here populate the T_mat_dis (lower triangle) with values of the
% residuals. NaN from original correction will be carried over and removed
T_mat_dis(get_indeces(N))=yresid;
T_mat_dis(find(isnan(T_mat_dis)==1))=0; % Here set negative edges to zero using the original T_mat!
%T_mat_dis(find(T_mat_dis<0))=0; % Here zero out the 
%T_mat_dis(find(D<Dthr))=0; % here also remove edges below min distance

[results_regress]=compute_SF(T_mat_dis,ind,1000);

% Here compute model fit if used polyfit
%SSresid=sum(yresid.^2);
%SStotal=(length(exp_vec)-1) * var(exp_vec);
%rsq=1 - SSresid/SStotal
%rsq_adj = 1 -SSresid/SStotal * (length(exp_vec)-1)/(length(exp_vec)-length(p))

%% Here determine the DS of the 136 genes using the DS_genes.csv file
ABA_dir='/nfs/zorba/ABA/Norm_March13/'
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/'
if ~exist('MA_resid')
    load([ABA_dir '/Science_Paper_MA_resid.mat'])
end
addpath([files_dir 'helper_functions'])

% Here grab the 136 gene names and extract from MA data
genes=importdata([files_dir '/136_genes.txt']);

% Here import the DS genes listed in Hawryzlyzc 2015
% and find the DS score for the 136 genes
DS_genes=importdata([files_dir '/DS_genes_local.csv'])
[c,ia,ib]=intersect(DS_genes.textdata(2:end,1),genes,'stable');
cons_DS=DS_genes.data(ia,4); % Index 4 is for cerebral cortex
genes=genes(ib);

%% Here run simulation for randomly selected clusters
%ABA_dir='/nfs/zorba/ABA/Norm_March13/'
ABA_dir='C:\Users\spiropan\Documents\Norm_March13\' % '/nfs/zorba/ABA/Norm_March13/'

if ~exist('MA_resid')
    load([ABA_dir 'Science_Paper_MA_resid.mat'])
end
addpath(['.\helper_functions'])

% Replicate primary analyses in Richiardi et. al. 
% Compute the total tissue similarity matrix, zero out negative edges and within-tissue edges
null_size=20
T_mat=corr(MA_resid); T_mat(find(T_mat<0))=0; 
T_mat(find(censor_mat==1))=0;
results=compute_SF(T_mat,ind,null_size);

% here create 100 random networks and replace the ind.W_all cell array and 
% recompute the SF each time
coords=load('Cortical_MNI_coords.csv'); 

clus_radius=[6:15];
for r = 1:length(clus_radius),
    for n=1:20
        ind_rand(n)=ind;
        [tmp]=random_clusters(coords,dat,clus_radius(r));
        %tmp=random_clusters_orig();
        ind_rand(n).W_all=tmp.W_all;
        ind_rand(n).Wi=tmp.Wi;
        [results_rand(n)]=compute_SF(T_mat,ind_rand(n),null_size);
        num_W_samples(n)=length(cat(2,tmp.W_all{[1:8,10:end]}));
        num_Wi_samples(n)=length(cat(2,tmp.W_all{[tmp.Wi]}));
    end
    
    % plot the distances for a sample simulated network as in Figure 1 of the reply from Richiardi (2017)
    [median_real(r),median_dist_rand(r),medians_dist_rand{r}]=plot_distances(ind,ind_rand,D);

    disp('percent of null networks with SF <0.05 uncorrected')
    for n=1:length(results_rand), pval_vec(n)=results_rand(n).pvalue; end
    perc_sig(r)=length(find(pval_vec<0.05))/length(pval_vec);
    median_W_samples(r)=median(num_W_samples);
    median_Wi_samples(r)=median(num_Wi_samples);
end

disp('correlation between percent significance and cluster size is:')
[r,p]=corr(perc_sig',clus_radius')

disp('correlation between clus radius and median number of W samples:')
[r,p]=corr(median_Wi_samples',clus_radius')

disp('correlation between percent significant SFs and median number of W samples:')
[r,p]=corr(median_W_samples',perc_sig')

disp('correlation between percent significant SFs and median number of Wi samples:')
[r,p]=corr(median_Wi_samples',perc_sig')

save(['Random_Results_nullsize_' int2str(null_size) '-' date],'results_rand','results')


