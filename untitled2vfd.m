clear all;close all;clc;

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
K_h = 1;  %Rician factors
K_g = 1;
finalResult3 =zeros(1, length(Power));
finalResult = zeros(1, length(Power));
finalresult2 = zeros(1,length(Power));
outPut=0;
 z=0.01;
for pIndex=1:30
    myz=0;
   
    for n=0:30
        
        for p=0:30
        
        myz = 4*exp(-(K_g+K_h))*K_h^(n)*K_g^(p)/(factorial(n)*factorial(p))^2 ...
            *z^(n+p+1)*besselk(n-p,2*z)+myz;
        
        end
        
    
    end
    z=z+5/30;
    finalResult(pIndex+1) = myz

end

plot(Power,finalResult)