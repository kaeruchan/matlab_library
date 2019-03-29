%initial
clear all;clc;
%parameter(Rayleigh)
m=1;

%b = @(k,gam) 2.*power(m,m)/(sqrt(2.*pi)*pow(gam,m)*gamma(m))* ...
% integral2(4.*pow(x,2.*n-1)*exp(-pow(y,2.)/2.-m*pow(x,2.)/gam-(2.*k+1.)*x*y-pow(x,2.)/2.),0,Inf,0,Inf);

%library
g_gam = @(gam) 0.5*(1+sqrt(1+2*m/gam)); %G(gamma)
b_1 = @(k,gam) sqrt(4.* gam/(2.* m + gam))* power(m,m) / ( power(gam, m)* ...
                                                  gamma(m))*1./(g_gam ...
                                                  + k); %B_k(m,1)
C_R =;



