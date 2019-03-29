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
F_gamma = zeros(1,length(P0));

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
        
  sum_main = 0;
        alpha = sqrt(2*V);
        beta = sqrt(2*(1+V)*p_Pj*(2^(2*R)*(1+x_k(n))-1)/(p_P0*Omega));
        c = 2*sqrt((1+V)*V*(M-1)/Omega);
        p = (1+V)/Omega;
        tilde_p = 2*p+beta^2;
        mu = M-1;
    
    %sum part
    for k=0:1:mu-1
        for j=0:1:mu-k-1
            sum_main = sum_main + (alpha/c)^(j+k) * (tilde_p/beta)^k * (-beta)^j ...
                * pochhammer(-mu+k+1 , j)/factorial(j)...
                * besseli(k+j,c*alpha*beta/tilde_p);
        end
    end
    
    
    % single channel
        F_gamma_d_k = 1 ...
                         - ((V+1)/Omega)^(M/2) * (V*(M-1))^(-(M-2)/2)*exp(-(M-1)*V) ...
                         * ((c/2)^(mu-1)*p^(-mu)*exp(c^2/(4*p))...
                         + 2/c*(c/tilde_p)^(mu)*exp((c^2/2-p*alpha^2)/tilde_p)...
                         * sum_main ...
                         - 1/p * (c/(2*p))^(mu-1)*exp(c^2/(4*p))...
                         * marcumq(beta*c/sqrt(2*p*tilde_p),alpha*sqrt(2*p/tilde_p),mu));
    % Totally
        F_gamma_d = 1 - (1- F_gamma_d_k)^M*(1- F_gamma_d_k);
 
     end
     F_gamma(Pindex) = F_gamma_d;
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