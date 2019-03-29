clean all;
close all;
clc;

 %from 0 to 90 degree

%rad_theta = 48*pi/180;% transform degree to rad
%r=0:10:3000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters library                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alpha beta gamma
alpha = 0.3; % suburban 0.1 Urban 0.3 Dense Urban 0.5
beta = 500; % suburban 750 Urban 500 Dense Urban 300
gamma = 15; % suburban 8 Urban 15 Dense Urban 20


C_a = [9.34e-01 2.30e-01 -2.25e-03 1.86e-05; ...,
       1.97e-02 2.44e-03 6.58e-06, 0; ...,
       -1.24e-04 -3.34e-06 0 0; ...
       2.73e-07 0 0 0];     
C_b = [1.17e-00 -7.56e-02 1.98e-03 -1.78e-05; ...
       -5.79e-03 1.81e-04 -1.65e-06 0; ...
       1.73e-05 -2.02e-07 0 0; ...
       -2.00e-08 0 0 0];
alpha_0 = 3.5;
alpha_pi = 2;
kappa_0 = 5;
kappa_pi = 15;
a1 = alpha_pi - alpha_0;
b1 = alpha_0;
sum_a2=0;
sum_b2=0;
for j=1:1:4
    for i=1:1:4
        sum_a2 = C_a(i,j) * (alpha*beta)^(i-1)*gamma^(j-1)+sum_a2;
    end
end
for j=1:1:4
    for i=1:1:4
        sum_b2 = C_b(i,j) * (alpha*beta)^(i-1)*gamma^(j-1)+sum_b2;
    end
end
a2 = sum_a2;
b2 = sum_b2;
a3 = kappa_0;
b3 = 2/pi*log(kappa_pi/kappa_0);



alpha_p = 0.5;
eta = 0.9;
%h = 500;
%r = 1;
R = 1;
h = 1:1:4000;
%Radious = h./tan(rad_theta);
kappa = 1.3;
Power= 20;
Pn = -150;
%m = 3;

p0 = 10.^(-3)*10.^(0.1*Power);
pn = 10.^(-3)*10.^(0.1*Pn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main program %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d=0;

%for rad_theta = 0:pi/180:pi/2;
    
    

r_op = zeros(1,length(h));
%for hlength=1:length(h)
for Radious =1:3000
    d_min = 9e+30;
    %for Radious=1:1:3000
    for rad_theta=pi/180:pi/180:pi*89/180
    %rad_theta = pi*73/180;
        d=0;
    for n=0:15
        for p=0:15
            %rad_theta = atan(h(hlength)/Radious);
            
            %rician factor
            V = a3 * exp(b3*rad_theta);
           
            V1 = a3 * b3 * exp(b3*rad_theta);
            

            %pathloss
            m = a1*1./(1+a2*exp(-b2*rad_theta))+b1;
            m1 = a1*a2*b2*exp(-b2*rad_theta)./(1+a2*exp(-b2*rad_theta)).^2;
            L= (2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0;
             epsilon = sqrt(L)*(csc(rad_theta))^m*(1+V) ...
                   *(Radious*tan(rad_theta))^m;   
            A = 4*(V+1)^2*exp(-2*V);
            A1 = (8*exp(-2*V)*(1+V)-8*exp(-2*V)*(1+V)^2)*V1;
            Z = (Radious*tan(rad_theta))^m*sqrt(L*(1+cot(rad_theta)^2)^m);
            Z1 = sqrt(L*(1+cot(rad_theta)^2)^m)...
               *(-(2*m*cot(rad_theta)*csc(rad_theta)^2)/(1+cot(rad_theta))+m1*log(1+cot(rad_theta)^2))...
               *(Radious*tan(rad_theta))^m/2 ...
               +sqrt(L*(1+cot(rad_theta)^2)^m)...
               *(m1*log(Radious*tan(rad_theta))+m*csc(rad_theta)*sec(rad_theta))*(Radious*tan(rad_theta))^m;
            G = V^(n+p)/(factorial(n)*factorial(p))^2*(V+1)^(n+p);
            G1 = (V1*V^(n+p-1)*(V+1)^(n+p)*(n+p)+V1*V^(n+p)*(V+1)^(n+p-1)*(n+p))...
                /(factorial(n)*factorial(p))^2;
            b = 2 * (V+1);
            M = MeijerG(-p,[],[n-p,0],-p-1,(b*Z)^2/4);
            d = -1/(factorial(n)*factorial(p))^2*2*V1*exp(-2*V)...
                *(V)^(n+p)*(1+V)^(2*p+2)*MeijerG(-p,[],[n-p,0],-1-p,epsilon^2)*epsilon^(2*(1+p))...
                +1/(factorial(n)*factorial(p))^2*V1*exp(-2*V)*V^(n+p-1)*(1+V)^(2*p+2)*(n+p)...
                *MeijerG(-p,[],[n-p,0],-1-p,epsilon^2)*epsilon^(2*(1+p)) ...
                +1/(factorial(n)*factorial(p))^2*2*V1*exp(-2*V)*V^(n+p)*(1+V)*(1+p)...
                *MeijerG(-p,[],[n-p,0],-1-p,epsilon^2)*epsilon^(2*(1+p))...
                +1/(factorial(n)*factorial(p))^2*exp(-2*V)*V^(n+p)*(1+V)*(2+2*p)...
                *MeijerG(-p,[],[n-p,0],-1-p,epsilon^2)*epsilon^(1+2*p)...
                *(L*(1+(cot(rad_theta))^2)^m*((-2*m*cot(rad_theta)*(csc(rad_theta))^2)/(1+(cot(rad_theta))^2)...
                  +m1*log(1+(cot(rad_theta))^2))*(Radious*tan(rad_theta))^m/(2*sqrt(L*(1+(cot(rad_theta))^2)^m))...
                  +sqrt(L*(1+(cot(rad_theta))^2)^m)*(m1*log(Radious*tan(rad_theta))+m*csc(rad_theta)*sec(rad_theta))*(Radious*tan(rad_theta))^m) ...
                +1/(L*(factorial(n)*factorial(p))^2)*exp(-2*V)*V^(n+p)*(1+V)^(-2)*(1+(cot(rad_theta))^2)^(-m)*(Radious*tan(rad_theta))^(-2*m) ...
                *epsilon^(2+2*p)*(2*V1*(1+V)*L*(1+(cot(rad_theta))^2)^m*(Radious*tan(rad_theta))^(2*m)...
                +(1+V)^2*L*(1+(cot(rad_theta))^2)^m*(-2*m*cot(rad_theta)*(csc(rad_theta))^2/(1+(cot(rad_theta))^2)+m1*log(1+(cot(rad_theta))^2))*(Radious*tan(rad_theta))^(2*m)...
                +(1+V)^2*L*(1+(cot(rad_theta))^2)^m*(2*m1*log(Radious*tan(rad_theta))+2*m*csc(rad_theta)*sec(rad_theta))*(Radious*tan(rad_theta))^(2*m))...
                *((-1-p)*MeijerG(-p,[],[n-p,0],[-1-p],epsilon^2)+2*besselk(n-p,2*epsilon)*epsilon^(n-p))+d;
                              
                
        end
    end
    d
    rad_theta
    if abs(d) < d_min
        d_min = abs(d)
        r_min = rad_theta;
    end        
    end
    r_op(Radious) = r_min
end
%end
%r_op
 %d_min
 %r_min

