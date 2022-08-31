clc;
clear;

%% Bladder
h=findobj(gca,'LineStyle',':','DisplayName','BLADDER');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x;y]);

f1=@(z) interp1(x,y,z);

x0=60; 
fprintf('BLADDER: \t V%.2f = %.2f\n', [x0;f1(x0)]);

%% Rectum
h=findobj(gca,'LineStyle',':','DisplayName','RECTUM');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x;y]);

f1=@(z) interp1(x,y,z);

x0=40; 
fprintf('RECTUM: \t V%.2f = %.2f\n', [x0;f1(x0)]);

%% CTV
h=findobj(gca,'LineStyle',':','DisplayName','CTV');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x1;y1]);

f1=@(z) interp1(x,y,z);

x0=78; 
fprintf('CTV: \t \t V%.2f = %.2f\n', [x0;f1(x0)]);

%% CTV
h=findobj(gca,'Type','line','DisplayName','CTV');
x1_tmp=h(2).XData;
y1_tmp=h(2).YData;

y1_first=find(y1_tmp==100,1,'last');
x1_tmp=x1_tmp(y1_first:end);
y1_tmp=y1_tmp(y1_first:end);
[y1,iy]=unique(y1_tmp);
x1=x1_tmp(iy);
f1=@(z) interp1(y1,x1,z);

x2_tmp=h(1).XData;
y2_tmp=h(1).YData;

%y2_first=find(y2_tmp==100,1,'last'); 
%x2_tmp=x2_tmp(y2_first+1:end);
%y2_tmp=y2_tmp(y2_first+1:end);
[y2,iy]=unique(y2_tmp);
x2=x2_tmp(iy);
f2=@(z) interp1(y2,x2,z);

x0=2; 
fprintf('CTV: \t \t V%.2f = [ %.2f - %.2f ]\n', [x0;f1(x0);f2(x0)]);

