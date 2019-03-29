clear all;
close all;
clc;

eta = 0.9;
h = 500;
R = 0.1:0.1:2;
Radious=100;
kappa = 1.3;
%r = Radious;
Power= 25;
Pn = -150;
m = 3;
K_h = 5;  %Rician factors
K_g = 5;

alpha3 = 0.2;
alpha5 = 0.4;
alpha8 = 0.6;

alphao = [0.825378,0.759787,0.712116,0.673708,0.641211 ...
           0.612903,0.587763,0.565125,0.544525,0.525627 ...
           0.508175,0.491972,0.476860,0.462711,0.449420 ...
           0.436898,0.425071,0.413876,0.403256,0.393165];

% finalResult = zeros(1, length(R));
% finalresult2 = zeros(1, length(R));
% finalresult3 = zeros(1, length(R));
% finalresult4 = zeros(1, length(R));
% finalresult5 = zeros(1, length(R));
% finalresult6 = zeros(1, length(R));
% finalresult8 = zeros(1, length(R));
% finalresult7 = zeros(1,length(R));
% for pIndex=1:length(Power)
%     outPut=0;
%     for r=0:1:(Radious)
%         myZ=0;
%         z1=0;
%         z2=0;
%         for n = 0 : 10
%             for p = 0 : 10
%                 
%                 %expressions
%                 %Power(pIndex)
%                 p0 = 10.^(-3)*10.^(0.1*Power(pIndex));
%                 pn = 10.^(-3)*10.^(0.1*Pn);
%                 L = (2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/(alpha)*pn/p0;
%                 Z = sqrt(L.*(h.^2+(r/(length(r))).^2).^m);
%                 a = 4*(K_h+1)*(K_g+1)*exp(-(K_g+K_h));
%                 G = K_g.^n*K_h.^p/(factorial(n)*factorial(p)).^2 ...
%                     *sqrt(((K_h+1)*(K_g+1)).^(n+p));
%                 b = 2*sqrt((K_h+1)*(K_g+1));
% 
% 
%                 z2 = a * G * Z^(2*p+2) *2^(n-p-2) *b^(p-n) ...
%                     * MeijerG(-p,[],[n-p,0],-p-1,b^2*Z^2/4)+z2;
% 
%             end
%         end
%         outPut = (r/(length(r))/ Radious^2) * z2/length(r)*4 + outPut
%     end
%     finalResult(pIndex) =  outPut
% end


%simulation

%for pIndex = 1:length(alpha)
for pIndex=1:length(R)
    %alphao = 1/(1+R(pIndex)*log(2));
    out=0;
    out2=0;
    out3=0;
    out4=0;
    out5=0;
    out6=0;
    out7=0;
    out8=0;
    for n=0:1000000      
        %expressions
       
        p0 = 10.^(-3)*10.^(0.1*Power);
        pn = 10.^(-3)*10.^(0.1*Pn);
        %g_k = ricernd(sqrt(K_g/(K_g+1)),sqrt(1/(K_g+1)/2));
        %h_k = ricernd(sqrt(K_h/(K_h+1)),sqrt(1/(K_g+1)/2));
        g_k = ricernd(sqrt(K_g/(K_g+1)),sqrt(1/(K_g+1)/2));
        h_k = ricernd(sqrt(K_h/(K_h+1)),sqrt(1/(K_h+1)/2));
        g_k_2 = ricernd(sqrt(K_g/(K_g+1)),sqrt(1/(K_g+1)/2));
        h_k_2 = ricernd(sqrt(K_h/(K_h+1)),sqrt(1/(K_g+1)/2));
        r_k = Radious;
        L3= (2^(R(pIndex)/(1-alpha3))-1)/kappa/eta*(1-alpha3)/alpha3*pn/p0;
        L5= (2^(R(pIndex)/(1-alpha5))-1)/kappa/eta*(1-alpha5)/alpha5*pn/p0;
        L8= (2^(R(pIndex)/(1-alpha8))-1)/kappa/eta*(1-alpha8)/alpha8*pn/p0;
        L0=(2^(R(pIndex)/(1-alphao(pIndex)))-1)/kappa/eta*(1-alphao(pIndex))/alphao(pIndex)*pn/p0;
        %a=(g_k*h_k)^2/(h^2+r_k^2)^m
        %b=(2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0
        if((g_k*h_k)^2/(h^2+r_k^2)^m < (L5))
            out =out+1;
        end
        if((g_k_2*h_k_2)^2/(h^2+r_k^2)^m < (L5))
            out2 =out2+1;
        end
        if((g_k_2*h_k_2)^2/(h^2+r_k^2)^m < (L3))
            out3 =out3+1;
        end
        if((g_k_2*h_k_2)^2/(h^2+r_k^2)^m < (L8))
            out4 =out4+1;
        end
        if((g_k_2*h_k_2)^2/(h^2+r_k^2)^m < (L8))
            out5 =out5+1;
        end
        if((g_k_2*h_k_2)^2/(h^2+r_k^2)^m < (L3))
            out6 =out6+1;
        end
        if((g_k_2*h_k_2)^2/(h^2+r_k^2)^m < (L0))
            out7 =out7+1;
        end
        if((g_k_2*h_k_2)^2/(h^2+r_k^2)^m < (L0))
            out8 =out8+1;
        end
    end
    finalResult(pIndex) = out/1000000
    finalresult2(pIndex) = out2/1000000
    finalresult3(pIndex) = out3/1000000
    finalresult4(pIndex) = out4/1000000
    finalresult5(pIndex) = out5/1000000
    finalresult6(pIndex) = out6/1000000
    finalresult7(pIndex) = out7/1000000
    finalresult8(pIndex) = out8/1000000
end
% %analysis alpha=0.5
% finalResult = [0.9923,0.9923,0.9923,0.9923,0.9923,0.9923,0.9922,0.9922,0.9920,0.9908, ...
%                0.9867,0.9756,0.9516,0.9085,0.8432,0.7569,0.6556,0.5481,0.4432,0.3477, ...
%                0.2657,0.1987,0.1461,0.1060,0.0762,0.0545,0.0389,0.0278,0.0199,0.0143, ...
%                0.0103];
% %simu. alpha =0.5
% finalresult2 = [1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,0.9999,0.9995,0.9975, ...
%                 0.9918,0.9775,0.9490,0.9024,0.8317,0.7438,0.6386,0.5307,0.4284,0.3339, ...
%                 0.2558,0.1929,0.1382,0.1016,0.0740,0.0512,0.0371,0.0262,0.0189,0.0130, ...
%                 0.0096];    
% 
% %simu alpha =0.3
% finalresult3 = [1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,0.9996, ...
%                 0.9982,0.9932,0.9812,0.9559,0.9137,0.8456,0.7614,0.6596,0.5533,0.4469, ...
%                 0.3515,0.2674,0.2018,0.1488,0.1090,0.0766,0.0549,0.0391,0.0294,0.0202, ...
%                 0.0139];
% 
% finalresult4 =  [1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,0.999970000000000,0.999820000000000,0.998830000000000,0.995250000000000,0.985970000000000,0.965500000000000,0.927300000000000,0.868890000000000,0.786410000000000,0.686920000000000,0.581740000000000,0.475410000000000,0.376640000000000,0.289380000000000,0.218050000000000,0.160700000000000,0.116700000000000,0.0847600000000000,0.0603300000000000,0.0448100000000000];
%             
% finalresult5 = [1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1.00001000000000,1,0.999970000000000,0.999620000000000,0.998870000000000,0.995040000000000,0.985880000000000,0.965120000000000,0.927460000000000,0.868450000000000,0.787810000000000,0.688310000000000,0.581040000000000,0.474450000000000,0.373400000000000,0.291280000000000,0.217600000000000,0.161370000000000,0.117450000000000,0.0854700000000000,0.0608100000000000,0.0439400000000000];            
%             
% %anal alpha =0.3
% finalresult6 = [0.9982,0.9982,0.9982,0.9982,0.9982,0.9982,0.9982,0.9982,0.9982,0.9980, ...
%                 0.9972,0.9938,0.9844,0.9632,0.9239,0.8625,0.7793,0.6797,0.5720,0.4653, ...
%                 0.3670,0.2818,0.2116,0.1560,0.1134,0.0817,0.0584,0.0417,0.0298,0.0213, ...
%                 0.0153];




clear title xlabel ylabel;
p1=semilogy(R,finalResult,'-');
ax = gca;
ax.FontSize=16;
ax.YLim = [0.0001,1];
grid on
p1.Color='Red';
p1.LineWidth=2;
xlabel('Transmission Rate','FontSize',16)
ylabel('Outage Probability','FontSize',16)


hold on
p2=semilogy(R,finalresult2,'v');
p2.Color='Red';
p2.MarkerSize=10;
%p2.LineWidth='2';

p3=semilogy(R,finalresult3,'*');
p3.Color='Blue';
p3.MarkerSize=10;
%p2.LineWidth='2';

p4=semilogy(R,finalresult4,'o');
p4.Color='Green';
p4.MarkerSize=10;
%p2.LineWidth='2';

p5=semilogy(R,finalresult5,'--');
p5.Color='Green';
p5.LineWidth=2;

p6=semilogy(R,finalresult6,':');
p6.Color='Blue';
p6.LineWidth=2;

p7=semilogy(R,finalresult7,'-.');
p7.Color='m';
p7.LineWidth=2;

p8=semilogy(R,finalresult8,'+');
p8.Color='m';
p8.MarkerSize=10;
hold off

lgd=legend([p6,p3,p1,p2,p5,p4,p7,p8],'Anal. \alpha=0.2','Simu. \alpha=0.2','Anal. \alpha=0.4', ...
    'Simu. \alpha=0.4','Anal. \alpha=0.6','Simu. \alpha=0.6','Anal. \alpha Optimized','Simu. \alpha Optimized');
lgd.FontSize=16;
lgd.Location='southeast';