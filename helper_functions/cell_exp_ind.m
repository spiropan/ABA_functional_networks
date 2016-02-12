function [out] = cell_exp_ind(input,expr,exactmatch)
% function [out] = cell_exp_ind(input,expr,exactmatch)
% Use this to return the row index (or column) of a
% cell array based on the expression being searched
% Will return all indeces whose string contains the 
% expr. Put exactmatch == 1 for exactmatch, 0 otherwise

out(1)=0;
count=1;

if nargin == 2
    exactmatch=0;
end

for h=1:length(input),
    if exactmatch == 0
        b=strfind(input{h},expr);
        %b=strfind(input(h,:),expr);
        %   if isempty(b{1}) == 0;
        if isempty(b) == 0;
            out(count)=h;
            count=count+1;
        end
    else
        b=strcmp(input{h},expr);
        if b==1
            out(count)=h;
            count=count+1;
        end
    end
end

