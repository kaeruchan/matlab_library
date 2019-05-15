close all;
clc;

% Parameter library
N_max = 100000;
M=5; %All relays
K=3; %Number of transmission
n=3; %Number of jammers
lambda=1; %Rayleigh Factor
Omega = 1; 
Pj = 5; %Power of Jammer (dBm)
P0 = 0:1:20; %Power of Signal (dBm)
gamma_k = 10; %threshold SNR (dB)
N0 = 0; %Power of noise (dBm)
q = 0.3; % Bernoulli Variance
kappa=0.3;

%initial
pr_anal_outage = zeros(1,length(P0));
pr_simu_outage1 = zeros(1,length(P0));
pr_simu_outage2 = zeros(1,length(P0));






for Pindex = 1 : length(P0)

    % Parameters (dbm --> watt)    
    p_Pj = 10^(Pj/10)*10^(-3);
    p_P0 = 10^(P0(Pindex)/10)*10^(-3);
    p_N0 = 10^(N0/10)*10^(-3);
    p_gamma_k = 10^(gamma_k/10);
    
    
    %initial
    R_sr = zeros(1,M);
    R_sr2 = zeros(1,M);
    R_sr2_2 = zeros(1,M-1);
    R_sr2_3 = zeros(1,M-2);
    R_sr_max = zeros(1,K);
    outage_counter = 0;
    outage_counter1 = 0;

    
    
    % Analysis
    
    
    % Simulation
    for N = 1:N_max
        
        %initial parameter
        gamma_main_simu1 = 0;
        %gamma_main_simu2 = 0;
        
        
        %Generate random number
        for Mnum = 1:M
            R_sr(Mnum) = (random('rayleigh',sqrt(1/(2*lambda))))^2;
            R_sr2(Mnum) = (random('rayleigh',sqrt(1/(2*lambda))))^2;
        end
        
        for Mnum = 1:M-1
            R_sr2_2(Mnum) = (random('rayleigh',sqrt(1/(2*lambda))))^2;
        end
        
         for Mnum = 1:M-2
            R_sr2_3(Mnum) = (random('rayleigh',sqrt(1/(2*lambda))))^2;
        end
        
        R_sr = sort(R_sr);
        R_sr2 = sort(R_sr2);
        
        
        for k = 1:K
            R_sr_max(k) = R_sr(M-k+1);
            gamma_main_simu1 = gamma_main_simu1 ...
                            + p_P0 ...
                            * min(R_sr_max(k),(random('rayleigh',sqrt(1/(2*lambda))))^2)...
                            /p_N0;
        end
        
        if (gamma_main_simu1 < p_gamma_k)
            outage_counter = outage_counter + 1;
        end
        
            gamma_main_simu2 =  p_P0 ...
                            * (min(max(R_sr2),random('rayleigh',sqrt(1/(2*lambda)))^2)...
                            + min(max(R_sr2_2),random('rayleigh',sqrt(1/(2*lambda)))^2)...
                            + min(max(R_sr2_3),random('rayleigh',sqrt(1/(2*lambda)))^2)) ...
                            /p_N0;
        
        if (gamma_main_simu2 < p_gamma_k)
            outage_counter1 = outage_counter1 + 1;
        end
        
    end
    pr_simu_outage1(Pindex) = outage_counter / N_max
    pr_simu_outage2(Pindex) = outage_counter1 / N_max

end


p1=semilogy(P0,pr_simu_outage1,'-');
ax = gca;
ax.FontSize=16;
ax.YLim = [0.001,1];
grid on
p1.Color = 'Red';
p1.LineWidth = 2;
xlabel('Transmission Power (P_0)','FontSize',16);
ylabel('Secrecy Outage Probability','FontSize',16);

hold on

p2 = semilogy(P0,pr_simu_outage2,'--v');
p2.Color = 'Red';
p2.MarkerSize = 10;
p2.LineWidth = 2;
