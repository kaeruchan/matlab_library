clear all;close all;clc;


Radious=100;
alpha = 0.5;
eta = 0.9;
h = 500;
%r = 100;
R = 1;
kappa = 1.3;
Power = 0:30;
Pn = -150;
m = 3;
K_h = 5;  %Rician factors
K_g = 5;

p=0;
n=0;


finalResult = zeros(1, length(Power));
finalresult2 = zeros(1,length(Power));

% Exponential implemented as Meijer G-function

;




for pIndex=1:length(Power)
    out=0;
    for n=0:100000      
        %expressions
       
        p0 = 10.^(-3)*10.^(0.1*Power(pIndex));
        pn = 10.^(-3)*10.^(0.1*Pn);
        g_k = ricernd(sqrt(K_g/(K_g+1)),sqrt(0.5/(K_g+1)));
        h_k = ricernd(sqrt(K_h/(K_h+1)),sqrt(0.5/(K_g+1)));
        %g_k = ricernd(sqrt(sqrt(K_g/(K_g+1))),sqrt(sqrt(1/(K_g+1)/2)));
        %h_k = ricernd(sqrt(K_h/(K_h+1)),sqrt(1/(K_h+1)/2));
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

plot(Power,finalresult2);