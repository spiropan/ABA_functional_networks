function [D]=euc_dis(pt1,pt2)
% function [D] = euc_dis(pt1,pt2)
    D=sqrt((pt1(1)-pt2(1))^2+(pt1(2)-pt2(2))^2+(pt1(3)-pt2(3))^2);
end