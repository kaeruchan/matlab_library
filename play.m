clear all;close all;clc;

x1 = 0:0.1:80;
x2 = 0:0.1:80;

y1 = ricepdf(x1/10,1,1);
y2 = ricepdf(x2/10,3,1);





plot (x1,y1);
%xlabel('”N—î','FontSize','12');
%ylabel('‰¿’l','FontSize','12');
hold on
plot (x2,y2);
hold off
