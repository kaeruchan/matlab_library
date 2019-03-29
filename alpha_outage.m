clear all;
close all;
clc;
%syms
%syms r;
%r = 0:5:100;
alpha = 0.05:0.05:0.95;
%alpha = 0.8;
eta = 0.9;
h = 500;
%r = 1;
R = 1;
Radious=100;
kappa = 1.3;
Power= 15;
Pn = -150;
m = 3;
K_h = 5;  %Rician factors
K_g = 5;

finalResult = zeros(1, length(alpha));
%finalresult2 = zeros(1,length(alpha));
%L = zeros(1,length(alpha));
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
%                 if abs(n-p)>1
%                     for k=0:abs(n-p)-1
%                         z1 = z1 + a/2*G*Z^(n+p+1)*(-1)^k*2^(-abs(n-p)-2*k) ...
%                             *(b*Z)^(2*k-abs(n-p))*factorial(abs(n-p)-k-1) ...
%                             /((2+2*k+n+p-abs(n-p))*factorial(k));
% %                 z2 = a * G * Z^(2*p+2) *2^(n-p-2) *b^(p-n) ...
% %                     * MeijerG(-p,[],[n-p,0],-p-1,b^2*Z^2/4)+z2;
%                     end
%                 end
%                 for k=0:10
%                     z2 = z2 + a*(-1)^n*G*b^(2*k+n)*Z^(2*n+2*k+p+1)*2^(n-2*k-1) ...
%                         *(b*Z)^(2*k+n) ...
%                         *(2+(2+2*k+2*n+p)*(psi(k+1)+psi(abs(n-p)+k+1)+2*log(2)-2*log(b*Z))) ...
%                         /(2+2*k+2*n+p)^2*factorial(k)*factorial(k+abs(n-p));                        
%                 end
%                 myZ = z1+z2+myZ;
%             end
%         end
%         outPut = (r/(length(r))/ Radious^2) * z2/length(r)*4 + outPut
%     end
%     finalResult(pIndex) =  outPut
% end


%simulation

for pIndex = 1:length(alpha)
%for pIndex=1:length(Power)
    out=0;
    out2=0;
    for n=0:100000      
        %expressions
       
        p0 = 10.^(-3)*10.^(0.1*Power);
        pn = 10.^(-3)*10.^(0.1*Pn);
        %g_k = ricernd(sqrt(K_g/(K_g+1)),sqrt(1/(K_g+1)/2));
        %h_k = ricernd(sqrt(K_h/(K_h+1)),sqrt(1/(K_g+1)/2));
        g_k = ricernd(sqrt(K_g/(K_g+1)),sqrt(1/(K_g+1)/2));
        h_k = ricernd(sqrt(K_h/(K_h+1)),sqrt(1/(K_h+1)/2));
        g_k_2 = ricernd(sqrt(K_g/(K_g+1)),sqrt(1/(K_g+1)/2));
        h_k_2 = ricernd(sqrt(K_h/(K_h+1)),sqrt(1/(K_g+1)/2));
        r_k = rand*Radious;
        L= (2.^(R./(1-alpha))-1)./kappa./eta.*(1-alpha)./alpha*pn./p0;
        %a=(g_k*h_k)^2/(h^2+r_k^2)^m
        %b=(2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0
        if(g_k*h_k)^2/(h^2+r_k^2)^m < (L(pIndex))
            out =out+1;
        end
        if(g_k_2*h_k_2)^2/(h^2+r_k^2)^m < (L(pIndex))
            out2 =out2+1;
        end
    end
    finalResult(pIndex) = out2/100000
    %finalresult2(pIndex) = out/100000
end

%power=20;
finalresult2 =[0.961070000000000,0.780240000000000,0.613540000000000,0.496580000000000,0.411290000000000,0.349960000000000,0.310160000000000,0.283890000000000,0.265850000000000,0.254860000000000,0.255360000000000,0.262760000000000,0.295200000000000,0.347170000000000,0.450950000000000,0.651550000000000,0.940990000000000,1.00001000000000,1.00001000000000];
finalresult3 =[0.568860000000000,0.287710000000000,0.175040000000000,0.123340000000000,0.0961900000000000,0.0759900000000000,0.0659900000000000,0.0592000000000000,0.0545100000000000,0.0524200000000000,0.0513800000000000,0.0536200000000000,0.0618900000000000,0.0784500000000000,0.108110000000000,0.199120000000000,0.501800000000000,0.997100000000000,1.00001000000000];
finalresult4 =[0.154950000000000,0.0606700000000000,0.0345900000000000,0.0232700000000000,0.0176000000000000,0.0152000000000000,0.0122200000000000,0.0106300000000000,0.0100400000000000,0.0103000000000000,0.0101700000000000,0.0104600000000000,0.0115500000000000,0.0134800000000000,0.0201200000000000,0.0376300000000000,0.126940000000000,0.814180000000000,1.00001000000000];

clear title xlabel ylabel;
p1=semilogy(alpha,finalResult,'-o');
ax = gca;
ax.FontSize=16;
ax.YLim = [0.001,1];
grid on
p1.Color='Red';
p1.DisplayName='Anal.';
p1.MarkerSize=10;
p1.LineWidth=2;
xlabel('\alpha','FontSize',16)
ylabel('Outage Probability','FontSize',16)


hold on
p2=semilogy(alpha,finalresult2,'--v');
p2.Color='m';
p2.MarkerSize=10;
p2.LineWidth=2;

p3=semilogy(alpha,finalresult3,':+');
p3.Color='Blue';
p3.MarkerSize=10;
p3.LineWidth=2;

p4=semilogy(alpha,finalresult4,'-.o');
p4.Color='Green';
p4.MarkerSize=10;
p4.LineWidth=2;

hold off

lgd=legend([p1,p2,p3,p4],'P_s = 15 dB','P_s = 20 dB','P_s = 25 dB','P_s = 30 dB');
lgd.FontSize=16;
lgd.Location='southwest';