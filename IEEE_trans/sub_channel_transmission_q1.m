%% Clear
clear all;close all;clc;

%% Parameter library
N_max = 10;
M=10; %All relays
K=7; %Number of transmission
n=3; %Number of jammers
Pj = 10; %Power of Jammer (dBm)
P0 = 25; %Power of Signal (dBm)
gamma_th = 10; %threshold SNR (dB)
N0 = 0; %Power of noise (dBm)
q = 0.1:0.05:0.9; % Bernoulli Variance
kappa=0.3;

% Channel -- Rician
V = 10; % Rician factor
Omega = 1; % Rician omega


%% Initial Parameters

% ORS
pr_anal_ORS_outage = zeros(1,length(q));
pr_simu_ORS_outage = zeros(1,length(q));

% ORSJ
pr_anal_ORSJ_outage = zeros(1,length(q));
pr_simu_ORSJ_outage = zeros(1,length(q));

% ORSMJ

pr_anal_ORSMJ_outage = zeros(1,length(q));
pr_simu_ORSMJ_outage = zeros(1,length(q));

% MRCMJ

pr_anal_MRCMJ_outage = zeros(1,length(q));
pr_simu_MRCMJ_outage = zeros(1,length(q));

% lines
pr_line_index = zeros(1,length(q));
% pr_line_MRCMJ = zeros(1,length(P0));
% pr_line_other = zeros(1,length(P0));
%% Program Start

for qindex = 1 : length(q)
    
    
    %% Parameters  (dBm --> Watt)
    p_Pj = 10^(Pj / 10)*10^(-3); % Jammer
    p_P0 = 10^(P0/10)*10^(-3); % Transmitter
    p_N0 = 10^(N0/10)*10^(-3); % Noise
    p_gamma_th = 10^(gamma_th/10); % threshold SNR
    
    %% lines
%     pr_line_MRCMJ(Pindex) = 1 - q;
%     pr_line_other(Pindex) = (1 - q)^2;
    
    %% Analysis
    
    % ORS
    
    
    % one hop outage
    MGF_ORS_one_hop = @(s) (1 - q(qindex) + exp(-V + V*(V+1)/(1+V+s*Omega))*q(qindex)*(1+V)/(1+V+s*Omega))^2 /s;
    
    
    
    pr_anal_ORS_outage(qindex) = euler_inversion(MGF_ORS_one_hop, p_gamma_th * p_N0/(kappa* p_P0));
    
    
    % ORSJ
    
    % sum part
    ORSJ_sum = 0;
    for k = 0:1:1
        for j = 0:1:(1-k)            
            ORSJ_sum =  ORSJ_sum + ((-1)^j * pochhammer(k-1,j))/(factorial(j)) ...
                       * (p_P0/(2*p_Pj*p_gamma_th))^((j+k)/2) ...
                       * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                       * besseli(j+k, 2*V*sqrt(2 * p_P0 * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th));
        end
    end
    pr_anal_ORSJ_outage(qindex) = 1 + 2 * (1-q(qindex)) * q(qindex) * ( ...  
                                  exp(-V*(p_P0 + n * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th)) ...  
                                  * p_Pj * p_gamma_th * (1 + V)/(p_P0 + p_Pj * p_gamma_th)...
                                  * besseli(0, 2*V*sqrt(p_P0 * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th))...
                                  - marcumq(sqrt(2*p_P0 * V/(p_P0 + p_Pj * p_gamma_th)), ...
                                            sqrt(2*V*p_Pj * p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 1)) ...
                                  + q(qindex)^2*(2^((n-1)/2)*exp(-V*(2*p_P0 + n * p_Pj * p_gamma_th)...
                                                         / (p_P0 + p_Pj * p_gamma_th)) ...
                                         * n * p_Pj^2 * p_gamma_th^2/(p_P0 + p_Pj * p_gamma_th)^2 ...
                                         * ORSJ_sum ...
                                         - marcumq(2*sqrt(p_P0 * V/(p_P0 + p_Pj * p_gamma_th)), ...
                                                   sqrt(2 * V * p_Pj * p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), ...
                                                   2));
              
    % ORSMJ
    
    %sum part 1
    ORSMJ_sum_1 = 0; %initial
    for k = 0:1:n-1
        for j = 0:1:n-k-1
            ORSMJ_sum_1 = ORSMJ_sum_1 ...
                        + (-1)^j*pochhammer(1+k-n,j)/factorial(j)...
                        * (p_P0*(V+1)/(p_Pj*p_gamma_th*Omega))^((j+k)/2) ...
                        * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                        * (n*Omega/(1+V))^((j+k+1)/2) ...
                        * besseli(1+j+k-n,2*V*sqrt(p_Pj * p_P0 * p_gamma_th * n) / (p_P0 + p_Pj * p_gamma_th));
        end
    end
    
    %sum part 2
    ORSMJ_sum_2 = 0; %initial
    for k = 0:1:n
        for j = 0:1:n-k
            ORSMJ_sum_2 = ORSMJ_sum_2 ...
                        + (-1)^j * pochhammer(k-n,j) * 2^((n-j-k-1)/2) / factorial(j) ...
                        * (p_P0 * (1+V)/(p_Pj * p_gamma_th * Omega))^((j+k)/2) ...
                        * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                        * (n * Omega / (V+1))^((j+k-1)/2) ...
                        * besseli(1+j+k-n,2 * V * sqrt(2 * p_Pj * p_P0 * p_gamma_th * n) / (p_P0 + p_Pj * p_gamma_th));
        end
    end
                               
    pr_anal_ORSMJ_outage(qindex) = 1 ...
                                 - (2 * (1-q(qindex)) * q(qindex) * marcumq( sqrt(2*p_P0*V/(p_P0+p_Pj*p_gamma_th)),sqrt(2*n*V*p_Pj*p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 1) ...
                                 - 2* (1-q(qindex)) * q(qindex) * (1+V)/sqrt(n * V) * exp(- V - (n-1) * p_Pj * V * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) ...
                                 * (p_P0* (1+V)/(p_Pj * p_gamma_th * Omega* n * V)) ^ ((n-1)/2) ...
                                 * (p_Pj * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) ^ n ...
                                 * (V * Omega / (1+V))^(n/2) ...
                                 * ORSMJ_sum_1 ...
                                 + q(qindex)^2 * marcumq(2 * sqrt(p_P0 * V / (p_P0 + p_Pj * p_gamma_th)), sqrt(2*n*V*p_Pj*p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 2) ...
                                 - q(qindex)^2 * p_Pj * p_gamma_th * sqrt(n/V) * exp(-2*V - (n-2)*p_Pj * V * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) / (p_P0 + p_Pj * p_gamma_th) ...
                                 * (p_P0 * (1+V) / (p_Pj * p_gamma_th * Omega * n * V)) ^ ((n-1)/2) ...
                                 * (p_Pj * p_gamma_th / (p_P0 + p_Pj * p_gamma_th))^n ...
                                 * (V * Omega / (1+V))^(n/2) ...
                                 * ORSMJ_sum_2);
                                    
                             
    % MRCMJ
    
    %sum part
    MRCMJ_sum = 0;
    for k = 0:1:n-1
        for j = 0:1:n-k-1
            MRCMJ_sum = MRCMJ_sum ...
                      + (-1)^j * pochhammer(1+k-n,j)/factorial(j) ...
                      * (p_P0 * (1+V)/(p_Pj * p_gamma_th * Omega))^((j+k)/2) ...
                      * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                      * (n * Omega /(1+V))^((j+k+1)/2) ...
                      * besseli(1+j+k-n, 2*V *sqrt(p_P0*p_Pj*n*p_gamma_th)/(p_P0 + p_Pj * p_gamma_th));
        end
    end
                     
    pr_anal_MRCMJ_outage(qindex) = 1 ...
                   - q(qindex)  ...
                   * marcumq(sqrt(2 * V*p_P0/(p_P0 + p_Pj * p_gamma_th)), sqrt(2 * n * V * p_Pj *p_gamma_th / (p_P0 + p_Pj *p_gamma_th)), 1) ...
                   + q(qindex)*(1+V)*exp(-(V + (n-1)*p_Pj*V*p_gamma_th/(p_P0 + p_Pj * p_gamma_th))) / (Omega * sqrt(n * V)) ...
                   * (p_P0 * (1+V)/(p_Pj * p_gamma_th * Omega * n * V))^((n-1)/2) ...
                   * (p_Pj * V * p_gamma_th/(p_P0 + p_Pj * p_gamma_th))^(n) ...
                   * (Omega/(V*(1+V)))^(n/2) ...
                   * MRCMJ_sum;
    
    
    %% Simulation 
    
    % ORS
    
    
    % Initial parameter
    outage_counter_ORS = 0;
    outage_counter_ORSJ = 0;
    outage_counter_ORSMJ = 0;
    outage_counter_MRCMJ = 0;
    
    
    for N = 1:N_max
        
        % ORS
        
        
        
        % Initial parameters
        
        % bernoulli distribution
        B1_ORS = random('Binomial',1,q(qindex));
        B2_ORS = random('Binomial',1,q(qindex));
        
        
        
        
        % f1
        f1_k_ORS = (random('rician', ...
                            sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
               
        % f2
        f2_k_ORS = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        gamma_inst_ORS = kappa * p_P0 * (B1_ORS * f1_k_ORS + B2_ORS)/p_N0;
        
        
        if gamma_inst_ORS < p_gamma_th
            outage_counter_ORS = outage_counter_ORS + 1;
        end
        
        
        % ORSJ
        % bernoulli distribution
        B1_ORSJ = random('Binomial',1,q(qindex));
        B2_ORSJ = random('Binomial',1,q(qindex));
        
        
        % f1
        f1_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        % f2
        f2_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        % v
        v_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        % gamma_inst
        gamma_inst_ORSJ = kappa * p_P0 * (B1_ORSJ * f1_ORSJ + B2_ORSJ * f2_ORSJ) / (kappa* p_Pj * v_ORSJ +p_N0);
        
        if gamma_inst_ORSJ < p_gamma_th
            outage_counter_ORSJ = outage_counter_ORSJ + 1;
        end   
        
        % ORSMJ
        
         % Initial parameters
         % bernoulli distribution
        B1_ORSMJ = random('Binomial',1,q(qindex));
        B2_ORSMJ = random('Binomial',1,q(qindex));
         
        
        % f1
        f1_ORSMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
     
        % f2
        f2_ORSMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        % v_n
        v_n_ORSMJ = 0;
        for i = 1:n
            v_n_ORSMJ = v_n_ORSMJ ...
                      + (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        gamma_inst_ORSMJ = kappa * p_P0 * (B1_ORSMJ * f1_ORSMJ + B2_ORSMJ * f2_ORSMJ) ...
                         / (p_N0 + kappa * p_Pj * v_n_ORSMJ);
        
        if gamma_inst_ORSMJ < p_gamma_th
            outage_counter_ORSMJ = outage_counter_ORSMJ + 1;
        end
        
        % MRCMJ
        
        % Initial parameters
        % bernoulli distribution
        B1_MRCMJ = random('Binomial',1,q(qindex));
        
        % f1
        f1_MRCMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        % v_n
        v_n_MRCMJ = 0;
        for i = 1:n
            v_n_MRCMJ = v_n_MRCMJ ...
                      + (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        gamma_inst_MRCMJ = kappa * p_P0 * B1_MRCMJ * f1_MRCMJ / (p_N0 + kappa * p_Pj * v_n_MRCMJ);
        
        if gamma_inst_MRCMJ < p_gamma_th
            outage_counter_MRCMJ = outage_counter_MRCMJ + 1;
        end
        
    end
    
    
    pr_simu_ORS_outage(qindex) = outage_counter_ORS / N_max;
    
    pr_simu_ORSJ_outage(qindex) = outage_counter_ORSJ /N_max;
    
    pr_simu_ORSMJ_outage(qindex) = outage_counter_ORSMJ / N_max;
    
    pr_simu_MRCMJ_outage(qindex) = outage_counter_MRCMJ / N_max
    
    
end


%% Plot

% ORS
p1=semilogy(q,pr_anal_ORS_outage,'-');
ax = gca;
ax.FontSize=16;
ax.YLim = [0.01,1];
grid on
p1.Color = 'Red';
p1.LineWidth = 2;
xlabel('Transmission Power (P_0)','FontSize',16);
ylabel('Eavesdropping Channel Transmission Outage Probability','FontSize',16);

hold on;

p2 = semilogy(q,pr_anal_ORS_outage,'v');
p2.MarkerSize = 10;
p2.Color = 'Red';

% ORSJ

p3=semilogy(q,pr_anal_ORSJ_outage,'-');
p3.LineWidth = 2;
p3.Color = 'Blue';

p4 = semilogy(q,pr_anal_ORSJ_outage,'^');
p4.MarkerSize = 10;
p4.Color = 'Blue';

% ORSMJ

p5=semilogy(q,pr_anal_ORSMJ_outage,'-');
p5.LineWidth = 2;
p5.Color = 'Cyan';

p6 = semilogy(q,pr_anal_ORSMJ_outage,'o');
p6.MarkerSize = 10;
p6.Color = 'Cyan';

% MRCMJ

p7=semilogy(q,pr_anal_MRCMJ_outage,'-');
p7.LineWidth = 2;
p7.Color = 'Magenta';

p8=semilogy(q,pr_anal_MRCMJ_outage,'+');
p8.MarkerSize = 10;
p8.Color = 'Magenta';

% p9=semilogy(q,pr_line_MRCMJ,'-.');
% p9.LineWidth = 2;
% p9.Color = 'Black';
% 
% p10=semilogy(q,pr_line_other,'--');
% p10.LineWidth = 2;
% p10.Color = 'Black';



% Line
p0=semilogy(q,pr_line_index,'-');
p0.LineWidth = 2;
p0.Color = 'Black';

% Brand
lgd=legend([p0,p2,p4,p6,p8],'Anal.','ORS','ORSJ','ORSMJ','MRCMJ');
lgd.Location = 'southwest';
lgd.FontSize = 14;

fname = '/users/shin/dropbox/programming/matlab';
saveas(p1,fullfile(fname,'sub_fig_j_q1'),'fig');