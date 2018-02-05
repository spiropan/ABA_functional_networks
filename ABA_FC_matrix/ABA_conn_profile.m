function [out] = ABA_conn_profile(input,varargin)
% function [out] = ABA_conn_profile(input, varargin)
% Use this to compute similarity in connectivity
% profiles in the input correlation or difference
% matrix. Basically between i and j this just removes 
% i and j from the vectors before computing correlation
% OPTIONS:
% pass 'weighted'   Pass matrix, same size as input, with weights to apply 
%                   to edges to estimate weighted correlations for
%                   connectivity profiles
% pass 'null'       Pass 1 to estimate null distribution and nonparametric p-vals
% pass 'dis_bins'   If passing null=1, pass dis_bins if want to permute
%                   edges within distance bins (that way it keeps
%                   dependence on distance in the null distribution
% pass 'dis_mat'     If passing dis_bins pass dis_mat same size as input          

% set defaults
dis_bins=[];
weighted=[];
dis_mat=[];
null=0;

% Parse vargin
i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'dis_bins'}
                dis_bins = varargin{i+1};
                i=i+2;
            case {'null'}
                null = varargin{i+1};
                i=i+2;
            case {'weighted'}
                weighted = varargin{i+1};
                i=i+2;
            case {'dis_mat'}
                dis_mat = varargin{i+1};
                i=i+2;
        end
    end
end


out=zeros(size(input,1),size(input,1));

% iterate over elements in the lower triangle of matrix
for j=1:size(input,2)-1
    for i=j+1:size(input,2)
        if i ~= j
            %disp(['working on regions ' int2str(i) ' and ' int2str(j)])
            tmp1=input(:,i);
            tmp2=input(:,j);
         
            tmp1([i,j])=[]; tmp2([i,j])=[]; % Remove self and connections between regionss
          
            % If doing weighted corrs
            if ~isempty(weighted)
                tmp1_weighted=weighted(:,i);
                tmp2_weighted=weighted(:,j);
                tmp1_weighted([i,j])=[]; tmp2_weighted([i,j])=[];
                  % Here take average for the weights and the inverse
                w=1./(tmp1_weighted+tmp2_weighted);
                % Compute weighted correlation
                temp=weightedcorrs([tmp1,tmp2],w);
                out(i,j)=temp(1,2);
                out(j,i)=temp(1,2);
            elseif null==1
                if ~isempty(dis_bins)
                    tmp2=tmp2(randperm(size(tmp2,1)));
                else
                    % Do the shuffling within distance bins
                    dis_vec=D(:,j); dis_vec([i,j])=[];
                    for d=1:length(dis_bins)-1   
                        in_ind=find(dis_vec > dis_bins(d) & dis_vec < dis_bins(d+1));
                        tmp2(in_ind)=tmp2(in_ind(randperm(length(in_ind))));
                    end
                end
                out(i,j)=corr(tmp1,tmp2);
                out(j,i)=out(i,j);                
            else
                % Compute unweighted
                out(i,j)=corr(tmp1,tmp2);
                out(j,i)=out(i,j);
            end
        end
    end
end


      
