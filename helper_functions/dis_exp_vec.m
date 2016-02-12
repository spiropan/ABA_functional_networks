function [out] = dis_exp_vec(T_mat,D_mat,censor_mat,ind)
% function [out] = dis_exp_vec(T_mat,D_mat,censor_mat,ind)
% Use this to create vectors of edges from tissue-tissue correlations and
% their corresponding distances
% out will have fields Wi_vec, TW_vec, DWi_vec, DTW_vec
% Pass the ind structure and has indeces for the Wi, T-W, edges etc.
% censor_mat is the matrix with ones for within-tissue class labels

N=size(T_mat,1);

T_mat(find(T_mat<0))=NaN; 
D_mat(find(T_mat<0))=NaN;

T_mat(find(censor_mat==1))=NaN;
D_mat(find(censor_mat==1))=NaN;

% here grab T_W indeces
disp('Working on T-W')
TW_mat=NaN(size(T_mat)); TW_mat(ind.W_all{9},ind.W_all{9})=T_mat(ind.W_all{9},ind.W_all{9});
TW_vec=TW_mat(get_indeces(N));
DTW_mat=D_mat; DTW_mat(find(isnan(TW_mat)==1))=NaN;
DTW_vec=DTW_mat(get_indeces(N));

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

% Here remove the NaN from the vectors
DWi_vec(find(isnan(DWi_vec)==1))=[];
Wi_vec(find(isnan(Wi_vec)==1))=[];
DTW_vec(find(isnan(DTW_vec)==1))=[];
TW_vec(find(isnan(TW_vec)==1))=[];

out.Wi_vec=Wi_vec;
out.DWi_vec=DWi_vec;
out.TW_vec=TW_vec;
out.DTW_vec=DTW_vec;
