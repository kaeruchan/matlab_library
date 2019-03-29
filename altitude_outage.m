clear all;
close all;
clc;
%syms
%syms r;
%r = 0:5:100;
alpha = 0.525627;
%alpha = 0.8;
eta = 0.9;
h = 0:50:1000;
%r = 1;
R = 1;
Radious=0:50:1000;
kappa = 1.3;
Power= 20;
Pn = -150;
m = 3;
K_h = 5;  %Rician factors
K_g = 5;

finalResult = zeros(length(Radious), length(h));

%simulation

% for rIndex = 1:length(Radious)
%     for hIndex = 1:length(h)
% %for pIndex=1:length(Power)
%         out=0;
%         out2=0;
%         for n=0:100000      
%         %expressions
%        
%             p0 = 10.^(-3)*10.^(0.1*Power);
%             pn = 10.^(-3)*10.^(0.1*Pn);
%         %g_k = ricernd(sqrt(K_g/(K_g+1)),sqrt(1/(K_g+1)/2));
%         %h_k = ricernd(sqrt(K_h/(K_h+1)),sqrt(1/(K_g+1)/2));
%             g_k = ricernd(sqrt(K_g/(K_g+1)),sqrt(1/(K_g+1)/2));
%             h_k = ricernd(sqrt(K_h/(K_h+1)),sqrt(1/(K_h+1)/2));
%             g_k_2 = ricernd(sqrt(K_g/(K_g+1)),sqrt(1/(K_g+1)/2));
%             h_k_2 = ricernd(sqrt(K_h/(K_h+1)),sqrt(1/(K_g+1)/2));
%             r_k = Radious(rIndex);
%             L= (2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0;
%         %a=(g_k*h_k)^2/(h^2+r_k^2)^m
%         %b=(2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0
%             if(g_k*h_k).^2./(h(hIndex).^2+r_k.^2).^m < (L)
%                out =out+1;
%             end
% %         if(g_k_2*h_k_2)^2/(h(hIndex)^2+r_k^2)^m < (L)
% %             out2 =out2+1;
% %         end
%         end
%         finalResult(rIndex,hIndex) = out/100000
%     %finalresult2(pIndex) = out/100000
%     end
% end
finalResult=[0,0,1.00000000000000e-05,1.00000000000000e-05,0.000190000000000000,0.000920000000000000,0.00347000000000000,0.0113300000000000,0.0363100000000000,0.103350000000000,0.241930000000000,0.458400000000000,0.698050000000000,0.879680000000000,0.969020000000000,0.995020000000000,0.999440000000000,1,1.00001000000000,1.00001000000000,1.00001000000000;0,0,0,0.000100000000000000,0.000300000000000000,0.00110000000000000,0.00382000000000000,0.0131200000000000,0.0398400000000000,0.109570000000000,0.250640000000000,0.470100000000000,0.707360000000000,0.883400000000000,0.970700000000000,0.995030000000000,0.999500000000000,0.999950000000000,1.00001000000000,1.00001000000000,1.00001000000000;0,1.00000000000000e-05,2.00000000000000e-05,0.000180000000000000,0.000460000000000000,0.00152000000000000,0.00481000000000000,0.0163900000000000,0.0471500000000000,0.123400000000000,0.279240000000000,0.501080000000000,0.734840000000000,0.898670000000000,0.974970000000000,0.996150000000000,0.999600000000000,0.999990000000000,1.00001000000000,1.00001000000000,1.00001000000000;1.00000000000000e-05,6.00000000000000e-05,0.000180000000000000,0.000380000000000000,0.000890000000000000,0.00257000000000000,0.00817000000000000,0.0238800000000000,0.0640800000000000,0.160960000000000,0.328030000000000,0.558750000000000,0.772350000000000,0.919530000000000,0.980780000000000,0.997220000000000,0.999750000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000;0.000260000000000000,0.000230000000000000,0.000330000000000000,0.000760000000000000,0.00225000000000000,0.00577000000000000,0.0148300000000000,0.0392800000000000,0.0978000000000000,0.213610000000000,0.402720000000000,0.631790000000000,0.828220000000000,0.943390000000000,0.987370000000000,0.998340000000000,0.999910000000000,1.0000000000000,1.0000000000000,1.00001000000000,1.0000000000000;0.00104000000000000,0.00117000000000000,0.00150000000000000,0.00274000000000000,0.00571000000000000,0.0127000000000000,0.0302200000000000,0.0707000000000000,0.152670000000000,0.300230000000000,0.500400000000000,0.717530000000000,0.879050000000000,0.964340000000000,0.993770000000000,0.999140000000000,0.999940000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;0.00338000000000000,0.00377000000000000,0.00520000000000000,0.00811000000000000,0.0148300000000000,0.0306400000000000,0.0603900000000000,0.125110000000000,0.239720000000000,0.413950000000000,0.621400000000000,0.805420000000000,0.927700000000000,0.981440000000000,0.996570000000000,0.999760000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;0.0114900000000000,0.0131100000000000,0.0164500000000000,0.0239600000000000,0.0396300000000000,0.0701000000000000,0.126420000000000,0.224340000000000,0.373280000000000,0.555510000000000,0.741390000000000,0.885590000000000,0.964280000000000,0.991260000000000,0.998960000000000,0.999950000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;0.0362300000000000,0.0397100000000000,0.0478500000000000,0.0646600000000000,0.0976900000000000,0.153240000000000,0.243970000000000,0.370960000000000,0.533350000000000,0.708560000000000,0.850750000000000,0.943030000000000,0.984060000000000,0.997200000000000,0.999730000000000,0.999990000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;0.104600000000000,0.108140000000000,0.127640000000000,0.160380000000000,0.217280000000000,0.296120000000000,0.415650000000000,0.555150000000000,0.709190000000000,0.839120000000000,0.931430000000000,0.978040000000000,0.994490000000000,0.999320000000000,0.999970000000000,1.000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;0.240570000000000,0.252190000000000,0.277890000000000,0.326510000000000,0.400480000000000,0.500870000000000,0.619100000000000,0.742860000000000,0.851670000000000,0.929820000000000,0.974920000000000,0.993260000000000,0.998700000000000,0.999890000000000,0.999990000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000;0.454930000000000,0.467740000000000,0.499790000000000,0.555980000000000,0.631830000000000,0.717000000000000,0.805510000000000,0.884930000000000,0.944060000000000,0.978270000000000,0.993560000000000,0.998650000000000,0.999780000000000,0.999990000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;0.696150000000000,0.708160000000000,0.733690000000000,0.775850000000000,0.825230000000000,0.880400000000000,0.926960000000000,0.963010000000000,0.984320000000000,0.994710000000000,0.998810000000000,0.999900000000000,1,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;0.881750000000000,0.885220000000000,0.898550000000000,0.918790000000000,0.944020000000000,0.964270000000000,0.980870000000000,0.991670000000000,0.997090000000000,0.999190000000000,0.999800000000000,0.999990000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.00001000000000,1.0000000000000;0.969150000000000,0.969560000000000,0.974690000000000,0.981440000000000,0.987710000000000,0.993440000000000,0.996580000000000,0.998830000000000,0.999590000000000,0.999950000000000,0.999970000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;0.995050000000000,0.994970000000000,0.996120000000000,0.997420000000000,0.998310000000000,0.999170000000000,0.999690000000000,0.999890000000000,0.999970000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.00001000000000,1.0000000000000,1.0000000000000;0.999540000000000,0.999590000000000,0.999670000000000,0.999750000000000,0.999920000000000,0.999970000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;0.999990000000000,0.999990000000000,0.999990000000000,1.0000000000000,1.0000000000000,1,1.0000000000000,1.00001000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.00001000000000,1.0000000000000,1.00001000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000;1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.000000000000,1.0000000000000,1.0000000000000,1.0000000000000;1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000,1.0000000000000];

%power=20;

clear title xlabel ylabel zlabel;
[h1,p1]=contour(h,Radious,finalResult,'ShowText','On');
%clabel(p1,x,'FontSize',15)
p1.TextStepMode = 'manual';
p1.TextStep=0.1;
p1.LineWidth = 2;
%p1.Fill='on';
ax = gca;
%clabel(p1,'FontSize',12);
p1.LevelStep = 0.1;
%p1.Fill ='on';
% ax.FontSize=16;
 ax.ZLim = [0,1];
 ax.XLim = [0,700];
 ax.YLim = [0,700];
grid on
% p1.Color='Red';
% p1.DisplayName='Anal.';
% p1.MarkerSize=10;
% p1.LineWidth=2;
ylabel('Radius (m)','FontSize',16)
xlabel('Altitude (m)','FontSize',16)
zlabel('Outage probability','FontSize',16)


% lgd=legend([p1,p2,p3],'r = 100 m','r = 200 m','r = 300 m');
% lgd.FontSize=16;
% lgd.Location='southeast';