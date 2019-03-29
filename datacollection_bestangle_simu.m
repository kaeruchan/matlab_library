%clean all; 
close all;
clc;

% angle_ana = [85,78,71,64,56,46,28,17,17,16 ...
%         16,16,16,15,15,15,15,15,15,15];
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % parameters library                                       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % alpha beta gamma
% alpha = 0.3; % suburban 0.1 Urban 0.3 Dense Urban 0.5
% beta = 500; % suburban 750 Urban 500 Dense Urban 300
% gamma = 15; % suburban 8 Urban 15 Dense Urban 20
% 
% 
% 
% C_a = [9.34e-01 2.30e-01 -2.25e-03 1.86e-05; ...,
%        1.97e-02 2.44e-03 6.58e-06, 0; ...,
%        -1.24e-04 -3.34e-06 0 0; ...
%        2.73e-07 0 0 0];     
% C_b = [1.17e-00 -7.56e-02 1.98e-03 -1.78e-05; ...
%        -5.79e-03 1.81e-04 -1.65e-06 0; ...
%        1.73e-05 -2.02e-07 0 0; ...
%        -2.00e-08 0 0 0];
% alpha_0 = 3.5;
% alpha_pi = 2;
% kappa_0 = 10^(0.1*5);
% kappa_pi = 10^(0.1*15);
% a1 = alpha_pi - alpha_0;
% b1 = alpha_0;
% sum_a2=0;
% sum_b2=0;
% for j=1:1:4
%     for i=1:1:4
%         sum_a2 = C_a(i,j) * (alpha*beta)^(i-1)*gamma^(j-1)+sum_a2;
%     end
% end
% for j=1:1:4
%     for i=1:1:4
%         sum_b2 = C_b(i,j) * (alpha*beta)^(i-1)*gamma^(j-1)+sum_b2;
%     end
% end
% a2 = sum_a2;
% b2 = sum_b2;
% a3 = kappa_0;
% b3 = 2/pi*log(kappa_pi/kappa_0);
% 
% 
% 
% alpha_p = 0.5;
% eta = 0.9;
% %h = 500;
% %r = 1;
% R = 1;
% %h = 1:1:4000;
% %Radious = h./tan(rad_theta);
% kappa = 1.3;
% Power= 20;
% Pn = -150;
% %m = 3;
% 
% p0 = 10.^(-3)*10.^(0.1*Power);
% pn = 10.^(-3)*10.^(0.1*Pn);    
%     
Radious = 0:10:410;
% r_op = zeros(1,length(Radious));
% 
% 
% 
%     
% for rIndex = 1:length(Radious)
%     outage_min = 1;
%     %for Radious=1:1:3000
%     for rad_theta=pi/180:pi/180:pi*89/180
%         counter=0;
%     %rad_theta = pi*73/180;
%         d=0;
%         for a=1:100000;
%             %rad_theta = atan(h(hlength)/Radious);
%             
%             %rician factor
%             V = a3 * exp(b3*rad_theta);
%            
%             V1 = a3 * b3 * exp(b3*rad_theta);
%             
% 
%             %pathloss
%             m = a1*1./(1+a2*exp(-b2*rad_theta))+b1;
%             m1 = a1*a2*b2*exp(-b2*rad_theta)./(1+a2*exp(-b2*rad_theta)).^2;
%             L= (2^(R/(1-alpha))-1)/kappa/eta*(1-alpha)/alpha*pn/p0;
%            
%             %channel fading coefficient
%             g = ricernd(sqrt(V/(V+1)),sqrt(1/(V+1)/2));
%             h = ricernd(sqrt(V/(V+1)),sqrt(1/(V+1)/2));
%             if (g*h) < sqrt(L*(Radious(rIndex)/cos(rad_theta))^(2*m))
%                    counter = counter+1;            
%             end
%         end        
%     outage = counter/100000;
%     if outage < outage_min
%         r_min = rad_theta*180/pi
%         outage_min = outage
%     end        
%     end
%     r_op(rIndex) = r_min
% end

r_op=[1,1,1,1,1,5,6,15,21,32,40,40,40,35,34,33,30,27,25,24,17,17,16,12,12,10,6,6,4,3,3,2,1,1,1,1,1,1,1,1,1,1];


p1=plot(Radious,r_op,'-o');
ax = gca;
ax.FontSize = 16;
ax.XLim = [0,400];
ax.YLim = [0,50];
grid on
p1.Color = 'Red';
p1.LineWidth =2;
xlabel('r_{max}[m]','FontSize',16);
ylabel('Optimum Angle [degree]','FontSize',16);
