% Clean all residual data
clear all;
close all;
clc;

% parameter library
N_max = 10000;
n_MAX = 10;
M=10; %The number of relay
V=5; %Rician Factor
Omega = 2; 
Pj = 5; %Power of Jammer (dBm)
P0 = 0:20; %Power of Signal (dBm)
gamma_k = 20; %threshold SNR (dBm)
N0 = 10; %Power of noise (dBm)
%beta = 0.5; %time ratio
q = 0.7; % Bernoulli Variance
R = 1; % secrecy rate

% initial parameter
pr_secrecy_analysis = zeros(1,length(P0));

for Pindex = 1:length(P0)
    
     % Parameters (dbm --> watt)    
    p_Pj = 10^(Pj/10)*10^(-3);
    p_P0 = 10^(P0(Pindex)/10)*10^(-3);
    p_gamma_k = 10^(gamma_k/10)*10^(-3);
    p_N0 = 10^(N0/10)*10^(-3);
    
    

    % Gauss-Laguerre library (15-point)

    x_k = [0.263560319718141 	    ... 1
           0.141340305910652e01   ... 2
           0.359642577104072e01   ... 3
           0.708581000585884e01   ... 4
           0.126408008442758e02   ... 5
          ];

    w_k = [0.521755610582809      ... 1
           0.398666811083176 	    ... 2
           0.759424496817076e-01  ... 3
           0.361175867992205e-02  ... 4
           0.233699723857762e-04  ... 5
           
          ];
   
      fx=0;
      pr=0;
     for n=1:5
        
         fx = fx ...
            + exp(x_k(n))*(2*exp(-((1+V)*x_k(n)/Omega+V))*(1-q)*q*(1+V)/Omega ...
            * besseli(0, 2*sqrt(V*(1+V)/Omega)*sqrt(x_k(n)))...
            + exp(-(2*V+(1+V)*x_k(n)/Omega))*q^2*sqrt((1+V)^3*x_k(n)/(2*V*Omega^3)) ...
            * besseli(1, 2*sqrt(2*V*(1+V)/Omega)*sqrt(x_k(n)))) ...
            * (1 - marcumq(sqrt(2*(M-1)*V),sqrt(2*(1+V)/Omega)*sqrt((p_P0*x_k(n))/(p_Pj*p_gamma_k)),M-1));
        
        pr = pr+fx*w_k(n);
     end
     pr_secrecy_analysis(Pindex) = 1-pr
end

%plot 
p1=semilogy(P0,pr_secrecy_analysis,'-');
ax = gca;
ax.FontSize=16;
ax.YLim = [0.05,1];
grid on
p1.Color='Black';
p1.LineWidth=2;
xlabel('Transmission Power (P_0)','FontSize',16)
ylabel('Sub-Channel Transmission Outage Probability','FontSize',16)
