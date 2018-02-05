function [ out ] = ABA_mantel_testBins( con_mat, exp_mat, dis_mat, varargin )
% function [ out ] = ABA_mantel_testBins( con_mat, exp_mat, dis_mat, vargin )
% Use this to run mantel test with permutation on connectivity and
% expression matrix
% OPTIONS:
% pass 'dis_bins'   (defaults is [0:16:144])
% pass 'con_thr'    (default is 0). Can be scalar or cell array matching
%                   number of distance bins
% pass 'null_size'  To estimate null distribution and nonparametric p-vals
% pass 'partial'    1 to adjust for distance (linear) 
% set defaults

dis_bins=[0:16:144];
con_thr=0;
null_size=0;
con_thr_cell=[];
partial = 0;
% Parse vargin
i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'dis_bins'}
                dis_bins = varargin{i+1};
                i=i+2;
            case {'con_thr'}
                if iscell(varargin{i+1})
                    con_thr_cell=varargin{i+1};
                else
                    con_thr_cell=[];
                    con_thr=varargin{i+1};
                end
                con_thr = varargin{i+1};
                i=i+2;
            case {'null_size'}
                null_size = varargin{i+1};
                i=i+2;
            case {'partial'}
                partial = varargin{i+1};
                i=i+2;
        end
    end
end

ind=get_indeces(size(con_mat,1));

full_con_vec=con_mat(ind);
exp_vec_orig=exp_mat(ind);
full_exp_vec=exp_vec_orig;
full_dis_vec=dis_mat(ind);

for i=1:null_size+1
    disp(['working on iteration ' int2str(i)])
    if i > 1           
        newind=randperm(size(exp_mat,1));        
        tmp=exp_mat(newind,newind);
        full_exp_vec=tmp(ind);
    end
    
    for d=1:length(dis_bins)-1      
        % here filter out the pairs below each min_distance
        
        % Apply a different threshold for each distance bin
        if ~isempty(con_thr_cell)
            con_thr=con_thr_cell{d};
        end
            
        if abs(con_thr) > 0
            if con_thr > 0
                in_ind=find(full_dis_vec>dis_bins(d) & full_dis_vec<dis_bins(d+1) & full_con_vec > con_thr);
            else
                in_ind=find(full_dis_vec>dis_bins(d) & full_dis_vec<dis_bins(d+1) & full_con_vec < con_thr);          
            end
        else 
           in_ind=find(full_dis_vec>dis_bins(d) & full_dis_vec<dis_bins(d+1));
        end
        
        dis_vec=full_dis_vec(in_ind); 
        exp_vec=full_exp_vec(in_ind); 
        con_vec=full_con_vec(in_ind);
             
        if i > 1
            exp_vec=exp_vec(randperm(length(exp_vec)));
        end
        
 
        %corrs(i,d)=corr(con_vec,exp_vec);
        
        try
            [corrs(i,d) p(i,d)]=corr(exp_vec,con_vec,'tail','both');
            % Here remove effects of distance from both connectivity and
            % expression
            %[B,dev,STATS_con]=glmfit(dis_vec,con_vec);
            %[B,dev,STATS_exp]=glmfit(dis_vec,exp_vec);
            if partial == 1
                [B,dev,STATS]=glmfit([dis_vec,exp_vec],con_vec);
                par_corrs(i,d)=B(2);
            end
            %[B,BINT,Res] = regress(dis_mat(ind),con_vec);
            %par_corrs(i,d)=corr(STATS_con.resid,STATS_exp.resid);
           
        catch
            corrs(i,d)=NaN;
            par_corrs(i,d)=NaN;
            p(i,d)=NaN;
        end
        %corrs(i)=corr(Res,exp_vec);
        
        % Only save these for the real values
        if i==1
            out.vals{i,d}.exp_vec=exp_vec;
            out.vals{i,d}.con_vec=con_vec;
            out.vals{i,d}.inds=in_ind;
        end
    end
end

for d=1:size(corrs,2)
    pval(1,d)=length(find(abs(corrs(2:end,d)) > abs(corrs(1,d))))/null_size;
    if partial ==1
        par_pval(1,d)=length(find(par_corrs(2:end,d) > par_corrs(1,d)))/null_size;
    else
        par_pval(1,d)=NaN;
    end
end

out.mant_R=corrs;
out.mant_pval=pval;
out.mant_parPval=p(1,:);
if partial == 1
    out.partial_mant_R=par_corrs;
    out.partial_mant_pval=par_pval;
end




