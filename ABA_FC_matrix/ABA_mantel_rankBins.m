function [] = ABA_mantel_rankBins( con_mat, AvgExp, dis_mat, name, varargin)
% function [] = ABA_mantel_rankBins( con_mat, AvgExp, dis_mat, name, con_thr, dis_bins)
% Use this for greedy search of genes that maximize ABA-FC Mantel correlation
% If ranking has already been done pass the name for the .mat
% file with the R variable (be sure to include .mat)
% con_thr if
% OPTIONS
% pass 'dis_bins' (default [32,144])
% pass 'con_thr'  (default is 0)
% set defaults

dis_bins=[32,144];
con_thr=0;
con_thr_cell=[];
varargin
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
                    con_thr_cell=varargin{i+1}
                else
                    con_thr_cell=[];
                    con_thr=varargin{i+1}
                end
                con_thr = varargin{i+1}
                i=i+2;
        end
    end
end


ind=get_indeces(size(con_mat,1));
num_bins=length(dis_bins)-1;

full_dis_vec=dis_mat(ind);
full_con_vec=con_mat(ind);

R=zeros(size(AvgExp,1),num_bins);
% Rows are the genes, go over each gene
if isempty(strfind(name,'.mat'))
    for d=1:num_bins,
        
        if ~isempty(con_thr_cell)
            con_thr=con_thr_cell{d}
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
        
        con_vec=full_con_vec(in_ind);
        dis_vec=full_dis_vec(in_ind);
        
        if isempty(con_vec)
            R(:,d)=NaN;
            continue
        end
        
        for i=1:size(AvgExp,1)         
            TmpExp=AvgExp;
            TmpExp(i,:)=[];
            TmpCorr=corr(TmpExp);
            
            exp_vec=TmpCorr(ind);
            exp_vec=exp_vec(in_ind);

            %[B,dev,STATS_con]=glmfit(dis_vec,con_vec);
            %[B,dev,STATS_exp]=glmfit(dis_vec,exp_vec);
            [B,dev,STATS]=glmfit([dis_vec,exp_vec],con_vec);
            
            %[B,BINT,Res] = regress(dis_mat(ind),con_vec);
            %par_corrs(i,d)=corr(STATS_con.resid,STATS_exp.resid);
            %par_corrs(i,d)=B(2);
            %[B,BINT,Res] = regress(dis_mat(ind),con_vec);
            
            R(i,d)=B(2);
            %R(i,d)=corr(STATS_con.resid,STATS_exp.resid); 
            %R(i,d)=corr(con_vec,exp_vec);
            
            %catch
                %R(i,d)=NaN;
            %end
            disp(['working on gene: ' int2str(i) '; distance bin: ' int2str(d) '; corr is ' num2str(R(i,d),'%0.15f')])
        end
    end
    save(name,'R','dis_bins')
else
    load(name)
    disp('Ok R mat successfully loaded')
    name(end-3:end)=[];
end

% Here rank the genes for each distance bin
if iscell(con_thr) 
    if con_thr{1} > 0 
        [Y,I]=sort(R,1,'ascend');
    elseif con_thr{1} < 0
        [Y,I]=sort(R,1,'descend');
    end
else
    [Y,I]=sort(R,1,'ascend');
end

size(Y)
size(I)
% Here plot the cumulative correlation for each distance bin?
Rcum=zeros(size(AvgExp,1),num_bins);
size(Rcum)
for d=1:num_bins,
    AvgExp_tmp=AvgExp(I(:,d),:);
    
    if abs(con_thr) > 0
        if con_thr > 0
            in_ind=find(full_dis_vec>dis_bins(d) & full_dis_vec<dis_bins(d+1) & full_con_vec > con_thr);
        else
            in_ind=find(full_dis_vec>dis_bins(d) & full_dis_vec<dis_bins(d+1) & full_con_vec < con_thr);
        end
    else
        in_ind=find(full_dis_vec>dis_bins(d) & full_dis_vec<dis_bins(d+1));
    end
    
    con_vec=full_con_vec(in_ind);
    dis_vec=full_dis_vec(in_ind);
    
    for i=1:size(AvgExp_tmp,1)       
        AvgExp_tmp(end,:)=[];
        try
            TmpCorr=corr(AvgExp_tmp);
            exp_vec=TmpCorr(ind);
            exp_vec=exp_vec(in_ind);
            
            % compute partial correlation
            [B,dev,STATS_con]=glmfit(dis_vec,con_vec);
            [B,dev,STATS_exp]=glmfit(dis_vec,exp_vec);
            
            Rcum(i,d)=corr(con_vec,exp_vec);
            disp('Ok able to compute cumulative Rcum')
        catch
            Rcum(i,d)=NaN;
        end
        disp(['working on cumulative gene: ' int2str(i) '; distance bin: ' int2str(d) '; corr is ' num2str(Rcum(i,d),'%0.15f')])
    end
end
save(name,'R','Rcum','Y','I','dis_bins')


