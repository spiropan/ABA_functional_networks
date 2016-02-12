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
    save([ABA_dir '/Science_Paper_MA_resid.mat'],'MA_resid','donor_ind','dat','gen','ABA_dir','files_dir')
else
    load([ABA_dir '/Science_Paper_MA_resid.mat'])
end

%% Extract the indeces for 13 networks, 4 networks of interest and from each subject
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA')
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
[results_resid]=compute_SF(T_mat,ind,1000);
 save([ABA_dir '/Science_Paper_MA_resid.mat'],'MA_resid','donor_ind','dat','gen','ABA_dir','files_dir',...
     'censor_mat','results_resid')
 
%% Here compute distances between all samples to compute distance matrix.
% columns 10-12 have MNI coordinates, 13 has network name
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA')
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
T_mat=corr(MA_resid);
T_mat(find(T_mat<0))=0; % Here zero out negative edges
T_mat(find(censor_mat==1))=0; % Here remove the tissue-tissue-edges too

for d=1:length(Dthr)
    NewT_mat=T_mat;
    NewT_mat(find(D<=Dthr(d)))=0; % Here zero out distances <X mm.
    [results_MinDis(d).out]=compute_SF(NewT_mat,ind,1000);
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
% out of network (T-W) edges
% Wi          = dark grey, 
% T-W         = light grey, 
% Order in ind.Wi is dDMN, Salience, Sensorimotor, Visuospatial
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA')
    load([ABA_dir '/Science_Paper_MA.mat']);
end
addpath([files_dir 'helper_functions'])

T_mat=corr(MA_resid); N=size(T_mat,1);

% Here zero out negative edges
T_mat(find(T_mat<0))=NaN; 
D_mat=D; D_mat(find(T_mat<0))=NaN;

% Here create a figure 
% here remove tissue-tissue connections
kill_figures
T_mat(find(censor_mat==1))=NaN;
D_mat(find(censor_mat==1))=NaN;

% here grab T_W indeces
disp('Working on T-W')
TW_mat=NaN(size(T_mat)); TW_mat(ind.W_all{9},ind.W_all{9})=T_mat(ind.W_all{9},ind.W_all{9});
TW_vec=TW_mat(get_indeces(N));
DTW_mat=D_mat; DTW_mat(find(isnan(TW_mat)==1))=NaN;
DTW_vec=DTW_mat(get_indeces(N));
%scatter(DTW_vec, TW_vec,'bo');
scatter(DTW_vec, TW_vec,[],[0.7,0.7,0.7],'filled');

% Here grab all the Wi indeces
disp('working on Wi')
Wi_ind=[];
for i=1:length(ind.Wi)
    Wi_ind=[Wi_ind ind.W_all{ind.Wi(i)}];
end
Wi_mat=NaN(size(T_mat)); Wi_mat(Wi_ind,Wi_ind)=T_mat(Wi_ind,Wi_ind);
Wi_vec=Wi_mat(get_indeces(N)); 
DWi_mat=D_mat; DWi_mat(find(isnan(Wi_mat)==1))=NaN;
DWi_vec=DWi_mat(get_indeces(N));

%hold on; scatter(DWi_vec,Wi_vec,'go');
hold on; scatter(DWi_vec, Wi_vec,[],[0.5,0.5,0.5],'filled');

% Here add the best fit line
DWi_vec(find(isnan(DWi_vec)==1))=[];
Wi_vec(find(isnan(Wi_vec)==1))=[];
DTW_vec(find(isnan(DTW_vec)==1))=[];
TW_vec(find(isnan(TW_vec)==1))=[];

tot_exp_vec=[Wi_vec; TW_vec]; tot_dis_vec=[DWi_vec; DTW_vec];
hold on; plot_linear(tot_dis_vec,tot_exp_vec,'k');
ylim([0,0.67])
xlim([0,180])
h_legend=legend('T-W','Wi');
set(h_legend,'fontsize',13);
xlabel('Euclidean Distance: mm','fontsize',16)
ylabel('Tissue-tissue correlation','fontsize',16)

[B,dev,stats]=glmfit(tot_dis_vec,tot_exp_vec);
[rho,pval]=corr(tot_dis_vec,tot_exp_vec);


%% Here can try a simple 2-sample ttest 
% This part will do 2-sample t-test 
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA')
    load([ABA_dir '/Science_Paper_MA.mat']);
end
addpath([files_dir 'helper_functions'])

% Here grab the 136 gene names and extract from MA data
genes=importdata([files_dir '/136_genes.txt']);
[c,ia,ib]=intersect(gen(:,3),genes);

Wi_ind=[];
for i=1:length(ind.Wi)
    Wi_ind=[Wi_ind ind.W_all{ind.Wi(i)}];
end
[p,h,ci,stats]=ttest2(MA_resid(:,Wi_ind)',MA_resid(:,ind.W_all{9})');
[h,p,ci,stats]=ttest2(MA_resid(:,ind.W_all{10})',MA_resid(:,ind.W_all{9})');
T_mat_small=corr(MA_resid(ia,:));
T_mat_small(find(T_mat_small<0))=0; % Here zero out negative edges
T_mat_small(find(censor_mat==1))=0;
%T_mat_small(find(D<20))=0;
[results_136genes]=compute_SF(T_mat_small,ind,1000);

%% here just determine distance-expression correlation using the 136 consensus genes
% 
exp_vec=T_mat_small(get_indeces(N));
D_mat=D; D_mat(find(T_mat_small<=0))=0;
dis_vec=D_mat(get_indeces(N));

exp_vec(find(exp_vec==0))=[];
dis_vec(find(dis_vec==0))=[];

%% Plot the effect of 136 consensus genes on mean tissue-tissue correlations binned by distance
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA')
    load([ABA_dir '/Science_Paper_MA.mat']);
end
addpath([files_dir 'helper_functions']);

% Here grab the 136 gene names and extract from MA data
genes=importdata([files_dir '/136_genes.txt']);
[c,ia,ib]=intersect(gen(:,3),genes);

T_mat_cons=corr(MA_resid(ia,:));
T_mat=corr(MA_resid);

consensus=dis_exp_vec(T_mat_cons,D,censor_mat,ind);
total=dis_exp_vec(T_mat,D,censor_mat,ind);

% Make the plots
% Make the distance break points
Dbr=[0:8:144]; Xtick={};
Wi_total=[]; TW_total=[]; Wi_cons=[]; TW_cons=[];
for i=1:length(Dbr)-1
    Wi_total(1,i)=median(total.Wi_vec(find(total.DWi_vec>Dbr(i) & total.DWi_vec<Dbr(i+1))));
    TW_total(1,i)=median(total.TW_vec(find(total.DTW_vec>Dbr(i) & total.DTW_vec<Dbr(i+1))));
    Wi_cons(1,i)=median(consensus.Wi_vec(find(consensus.DWi_vec>Dbr(i) & consensus.DWi_vec<Dbr(i+1))));
    TW_cons(1,i)=median(consensus.TW_vec(find(consensus.DTW_vec>Dbr(i) & consensus.DTW_vec<Dbr(i+1))));
    %Xtick=[Xtick [int2str(Dbr(i)) ':' int2str(Dbr(i+1))]];
end
kill_figures
%plot(Wi_total,'ko-')
%hold on; plot(TW_total,'k*-')
%hold on; plot(Wi_cons,'ks-')
%hold on; plot(TW_cons,'k.-')

plot(Wi_total,'o-','Color',[0,0,0])
hold on; plot(TW_total,'*-','Color',[0,0,0,])
hold on; plot(Wi_cons,'o-','Color',[0.5,0.5,0.5])
hold on; plot(TW_cons,'*-','Color',[0.5,0.5,0.5])

h_legend=legend('Wi: all genes','T-W: all genes','Wi: consensus','T-W: consensus')
set(h_legend,'fontsize',12);
xlabel('Distance Bins: [0:8:144] mm','fontsize',16)
ylabel('Median tissue-tissue correlation','fontsize',16)
xlim([0,19])
ylim([0.05,0.16])
% Here are the edges that show the spike in the graph (128-136)
% in Wi for consensus genes
dmin=128; dmax=136;
%length(consensus.Wi_vec(find(consensus.DWi_vec > dmin & consensus.DWi_vec < dmax)))
%std(consensus.Wi_vec(find(consensus.DWi_vec > dmin & consensus.DWi_vec < dmax)))
%plot(consensus.Wi_vec(find(consensus.DWi_vec > dmin & consensus.DWi_vec < dmax)))

%% Here plot strength fraction as a function of (binned) distances
ABA_dir='/nfs/zorba/ABA/Norm_March13/';
files_dir='/home/spantazatos/Dropbox/Postdoc/Science_Commentary/';
if ~exist('MA')
    load([ABA_dir '/Science_Paper_MA.mat']);
end
addpath([files_dir 'helper_functions'])

% Standard preprocessing below, remove negative edges and within-tissue
% edges
T_mat=corr(MA_resid);
T_mat(find(T_mat<0))=0;
T_mat(find(censor_mat==1))=0;

% Make the plots
% Make the distance break points
Dbr=[0:8:144]; Xtick={}; null_size=1000;
SF=zeros(null_size+1,length(Dbr)-1);
for i=1:length(Dbr)-1
    tmp=zeros(size(T_mat));
    tmp(find(D>Dbr(i) & D<Dbr(i+1)))=T_mat(find(D>Dbr(i) & D<Dbr(i+1)));
    out(i)=compute_SF(tmp, ind, null_size);
    SF(1,i)=out(i).real_SF;
    SF(2:end,i)=out(i).null_SF';
    %Xtick=[Xtick [int2str(Dbr(i)) ':' int2str(Dbr(i+1))]];
end
results.accuracy=SF;
kill_figures
save('to_plot','results');
spider_SVM_plot({'to_plot.mat'},[1:length(Dbr)-1],{'Real SF','Null_SF'})
h_legend=legend('Real SF','Null SF')
set(h_legend,'fontsize',14)
%h=gca;
%set(h,'XTickLabel',{'','Tissue','<4 mm','<8 mm','<12 mm','<16 mm','<20 mm','<24 mm'},'fontsize',12)
ylabel('Strength Fraction','fontsize',14)
xlabel('Distance: [0:8:144] mm bins','fontsize',14)

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
DS_genes=importdata([files_dir '/DS_genes.csv'])
[c,ia,ib]=intersect(DS_genes.textdata(2:end,1),genes,'stable');
cons_DS=DS_genes.data(ia,1);
genes=genes(ib);



