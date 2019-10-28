function [out] = compute_SF(T_mat, ind, null_size)
% function [out, perm_ind] = compute_SF(T_mat, ind, null_size)
% This will compute strength fraction as defined in Richiardi et. al
% ind is a structure with:
%    ind.W_all{n} = indeces for 13 networks plus Z_restOfBrain
%    ind.Wi= indeces of ind.W_all for dDMN, Salience, Sensorimotor, Visuospatial
%    ind.W = indeces of ind.W_all for 13 networks (i.e. just excluding
%    Z_restof Brain

real_SF = SF_calc(T_mat,ind);

for i=1:null_size
    %disp(['working on null iteration: ' int2str(i)]);
    perm_ind=[];
    % here permute the samples columns of MA keeping structure of brain
    for d=1:length(ind.donor_ind)
        tmp=ind.donor_ind{d}(randperm(length(ind.donor_ind{d})));
        perm_ind=[perm_ind tmp];
    end
    
    Tmat_perm=T_mat(perm_ind,perm_ind);
    null_SF(i)=SF_calc(Tmat_perm,ind);
end

pvalue=length(find(null_SF > real_SF))/null_size;

out.real_SF=real_SF;
out.null_SF=null_SF;
out.pvalue=pvalue;
end

function [SF] = SF_calc(T_mat,ind)
% Will use helper function get_indeces(N) to get lower traingle indeces fo
% NxN correlation matrix

% First compute T over the whole graph
N=size(T_mat,1); T=sum(T_mat(get_indeces(N)));
% Here compute the metrics W
W=0; % initialize W
for i=1:length(ind.W)
    % ind.W_all{ind.W(i)} will have sample indeces for each network
    W_mat=T_mat(ind.W_all{ind.W(i)},ind.W_all{ind.W(i)}); N=size(W_mat,1); W=W + sum(W_mat(get_indeces(N)));
end
% Here compute the netric Wi
Wi=0; %initialize Wi
for i=1:length(ind.Wi)
    % ind.W_all{ind.Wi(i)} will have sample indeces for each of four
    % networks
    Wi_mat=T_mat(ind.W_all{ind.Wi(i)},ind.W_all{ind.Wi(i)}); N=size(Wi_mat,1); Wi=Wi + sum(Wi_mat(get_indeces(N)));
end
SF=Wi/(T-W);
end