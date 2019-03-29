clear all; close all; clc;
alpha = 0.5;
eta = 0.9;
h = 500;
%r = 1;
R = 1;
Radious=100;
kappa = 1.3;
Power= 0:30;
Pn = -150;
m = 3;
K_h = 5;  %Rician factors
K_g = 5;

finalresult2 = zeros(1,length(Power));

out=0;


for pIndex=1:length(Power)
    out=0;
    for n=0:100000      
        %expressions
       
        p0 = 10.^(-3)*10.^(0.1*Power(pIndex));
        pn = 10.^(-3)*10.^(0.1*Pn);
        g_k = ricernd(sqrt(K_g/(K_g+1)),sqrt(0.5/(K_g+1)));
        h_k = ricernd(sqrt(K_h/(K_h+1)),sqrt(0.5/(K_g+1)));
        r_k = rand*Radious;
        L= (2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0;
        %a=(g_k*h_k)^2/(h^2+r_k^2)^m
        %b=(2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0
        if((g_k*h_k)^2/(h^2+r_k^2)^m < (L))
            out =out+1;
        end
    end
    finalresult2(pIndex) = out/100000
end

semilogy(Power,finalresult2)