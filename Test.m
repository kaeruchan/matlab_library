clear all;close all;clc;
%syms
%syms r;
%r = 0:5:100;

alpha = 0.5;
eta = 0.9;
h = 100;
%r = 1;
R = 1;
Radious=50;
kappa = 1.3;
Power= 0:30;
Pn = -200;
m = 3;
K_h = 5;  %Rician factors
K_g = 5;

finalResult = zeros(1, length(Power));
% finalresult2 = zeros(1, length(Power));
outPut= 0;
for pIndex=1:length(Power)
    for r=0:2:Radious
        myZ=0;
        for n = 0 : 10
            for p = 0 : 10
                
                %expressions
                p0 = 10.^(-3)*10.^(0.1*Power(pIndex));
                pn = 10.^(-3)*10.^(0.1*Pn);
                L = (2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/(alpha)*pn/p0;
                Z = sqrt(L.*(h.^2+r.^2).^m);
                a = 4*(K_h+1)*(K_g+1)*exp(-(K_g+K_h));
                G = K_g*K_h/(factorial(n)*factorial(p)).^2*sqrt(((K_h+1)*(K_g+1)).^(n+p));
                b = 2*sqrt((K_h+1)*(K_g+1));
                
                for z = 0:2:Z
                   
                 %z=0;
                 myZ = a.*G.*(z)^(n+p+1)*besselk(n-p,b*z)+ myZ;
                end
            end
        end
        outPut = (r/ Radious^2) * myZ + outPut;
    end
    finalResult(pIndex) = outPut
    
%     %simulation
%     out=0;
%     for n=0:100000      
%         %expressions
%         p0 = 10.^3*10.^(0.1*Power(pIndex));
%         pn = 10.^3*10.^(0.1*Pn);
%         g_k = ricernd(K_g,1);
%         h_k = ricernd(K_h,1);
%         r_k = rand*Radious;
%         %a=(g_k*h_k)^2/(h^2+r_k^2)^m
%         %b=(2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0
%         if((g_k*h_k)^2/(h^2+r_k^2)^m < (2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0)
%             out =out+1;
%         end
%     end
%     finalresult2(pIndex) = out/100000
end




semilogy(Power,finalResult)
hold on
semilogy(Power,finalresult2)
hold off
