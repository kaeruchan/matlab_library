clear all;close all;clc;
%syms
%syms r;
%r = 0:5:100;

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
finalResult3 =zeros(1, length(Power));
finalResult = zeros(1, length(Power));
finalresult2 = zeros(1,length(Power));
outPut=0;
for pIndex=1:length(Power)
    z3=0;
    for r=0:2:Radious
        z2=0;
        for n = 0 : 30
            for p = 0 : 30
                
                %expressions
                p0 = 10.^3*10.^(0.1*Power(pIndex));
                pn = 10.^3*10.^(0.1*Pn);
                L = (2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/(alpha)*pn/p0;
                Z = sqrt(L.*(h.^2+r.^2).^m);
                a = 4*(K_h+1)*(K_g+1)*exp(-(K_g+K_h));
                G = K_g^(2*n)*K_h^(2*p)/(factorial(n)*factorial(p)).^2*sqrt(((K_h+1)*(K_g+1)).^(n+p));
                b = 2*sqrt((K_h+1)*(K_g^2+1));
                
               
                if abs(n-p) >1
                    z11=0;
                    for k1 = 0:abs(n-p)-1
                        z11 = a/2*G*Z^(n+p+1)*(-1)^k1*2^(-abs(n-p)+2*k1)*(b*Z)^(2*k1-abs(n-p)) ...
                        *factorial(abs(n-p)-k1-1)/((2+2*k1+n+p-abs(n-p))*factorial(k1))+z11;
                    end
                else z11=0;
                end
                    z12=0;
                    for k2 = 0:10
                        z12 = (-1)^n * a * G *b^(2*k2+n)*Z^(2*n+2*k2+p+1)*2^(n-2*k2-1)*(b*Z)^(2*k2+n) ...
                            *(2+(2+2*k2+2*n+p)*(psi(k2+1)+psi(abs(n-p)+k2+1)+2*log(2)-2*log(b*Z))) ...
                            /((2 + 2*k2 + 2*n +p)^2*factorial(k2)*factorial(k2+abs(n-p)))+z12;
                    end
                    z2 = z11+ z12+z2;
              end
              
            
        end
        outPut = (r/ Radious^2) * z2 + outPut;
        %outPut = 1 / Radious * z2 + outPut;
    end
    %finalResult3(pIndex) = 1 - 2* sqrt((h^2+Radious^2)^m * L)*besselk(1,2*sqrt((h^2+Radious^2)^m*L));
    finalResult(pIndex) =  outPut
end


%simulation


for pIndex=1:length(Power)
    out=0;
    for n=0:100000      
        %expressions
       
        p0 = 10.^3*10.^(0.1*Power(pIndex));
        pn = 10.^3*10.^(0.1*Pn);
        g_k = ricernd(sqrt((K_g)/((K_g+1))),sqrt(1/(K_g+1)/(2)));
        h_k = ricernd(sqrt((K_h)/((K_h+1))),sqrt(1/(K_g+1)/(2)));
        r_k = Radious;
        
        L= (2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0;
        %a=(g_k*h_k)^2/(h^2+r_k^2)^m
        %b=(2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0
        if  (g_k*h_k)^2/(h^2+r_k^2)^m < (L)
            out =out+1;
        end
    end
    finalresult2(pIndex) = out/100000
end


plot(Power,finalResult3)
hold on
plot(Power,finalresult2)
hold off