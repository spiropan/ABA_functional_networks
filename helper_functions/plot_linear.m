function plot_linear(x,y,col)
% usage: plot_linear(X,Y,'b')
range=abs(max(x)-min(x))

padding=0.10*range

xFitting=[min(x)-padding:0.001:max(x)+padding]

[coefs]=polyfit(x, y, 1);
yFitted = polyval(coefs, xFitting);
colstr = sprintf('o%s',col);
%plot(x,y,colstr,x,i2,col,'MarkerSize',7);

% uncomment below if want to plot all points
%plot(xFitting,yFitted,col,x,y,'k*','MarkerSize',5);
plot(xFitting,yFitted,col,'LineWidth',3);
xlim([min(x)-padding,max(x)+padding]);
