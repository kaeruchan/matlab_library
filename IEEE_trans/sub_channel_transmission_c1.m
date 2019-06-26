%% Clear
clear all;close all;clc;

%% Parameter library
N_max = 100000;
M=10; %All relays
K=7; %Number of transmission
n=3; %Number of jammers
Pj = 10; %Power of Jammer (dBm)
P0 = 0:2:40; %Power of Signal (dBm)
gamma_th = 10; %threshold SNR (dB)
N0 = 0; %Power of noise (dBm)


q = 0.3; % Bernoulli Variance
q2 = 0.5;
q3 = 0.7;


kappa=0.3;

% Channel -- Rician
V = 10; % Rician factor
Omega = 1; % Rician omega


%% Initial Parameters

% ORS
pr_anal_ORS_outage = zeros(1,length(P0));
pr_simu_ORS_outage = zeros(1,length(P0));
pr_anal_ORS2_outage = zeros(1,length(P0));
pr_simu_ORS2_outage = zeros(1,length(P0));
pr_anal_ORS3_outage = zeros(1,length(P0));
pr_simu_ORS3_outage = zeros(1,length(P0));

% ORSJ
pr_anal_ORSJ_outage = zeros(1,length(P0));
pr_simu_ORSJ_outage = zeros(1,length(P0));
pr_anal_ORSJ2_outage = zeros(1,length(P0));
pr_simu_ORSJ2_outage = zeros(1,length(P0));
pr_anal_ORSJ3_outage = zeros(1,length(P0));
pr_simu_ORSJ3_outage = zeros(1,length(P0));

% ORSMJ

pr_anal_ORSMJ_outage = zeros(1,length(P0));
pr_simu_ORSMJ_outage = zeros(1,length(P0));
pr_anal_ORSMJ2_outage = zeros(1,length(P0));
pr_simu_ORSMJ2_outage = zeros(1,length(P0));
pr_anal_ORSMJ3_outage = zeros(1,length(P0));
pr_simu_ORSMJ3_outage = zeros(1,length(P0));


% MRCMJ

pr_anal_MRCMJ_outage = zeros(1,length(P0));
pr_simu_MRCMJ_outage = zeros(1,length(P0));
pr_anal_MRCMJ2_outage = zeros(1,length(P0));
pr_simu_MRCMJ2_outage = zeros(1,length(P0));
pr_anal_MRCMJ3_outage = zeros(1,length(P0));
pr_simu_MRCMJ3_outage = zeros(1,length(P0));


% lines
pr_line_index = zeros(1,length(P0));
pr_line_index2 = zeros(1,length(P0));
pr_line_index3 = zeros(1,length(P0));
pr_line_MRCMJ = zeros(1,length(P0));
pr_line_other = zeros(1,length(P0));
pr_line_MRCMJ2 = zeros(1,length(P0));
pr_line_other2 = zeros(1,length(P0));
pr_line_MRCMJ3 = zeros(1,length(P0));
pr_line_other3 = zeros(1,length(P0));

%% Program Start

for Pindex = 1 : length(P0)
    
    
    %% Parameters  (dBm --> Watt)
    p_Pj = 10^(Pj / 10)*10^(-3); % Jammer
    p_P0 = 10^(P0(Pindex)/10)*10^(-3); % Transmitter
    p_N0 = 10^(N0/10)*10^(-3); % Noise
    p_gamma_th = 10^(gamma_th/10); % threshold SNR
    
    %% lines
    pr_line_MRCMJ(Pindex) = 1 - q;
    pr_line_other(Pindex) = (1 - q)^2;
    pr_line_MRCMJ2(Pindex) = 1 - q2;
    pr_line_other2(Pindex) = (1 - q2)^2;
    pr_line_MRCMJ3(Pindex) = 1 - q3;
    pr_line_other3(Pindex) = (1 - q3)^2;
    %% Analysis
    
    % ORS
    
    
    % one hop outage
    MGF_ORS_one_hop = @(s) (1 - q + exp(-V + V*(V+1)/(1+V+s*Omega))*q*(1+V)/(1+V+s*Omega))^2 /s;
    MGF_ORS2_one_hop = @(s) (1 - q2 + exp(-V + V*(V+1)/(1+V+s*Omega))*q2*(1+V)/(1+V+s*Omega))^2 /s;
    MGF_ORS3_one_hop = @(s) (1 - q3 + exp(-V + V*(V+1)/(1+V+s*Omega))*q3*(1+V)/(1+V+s*Omega))^2 /s;
    
    pr_anal_ORS_outage(Pindex) = euler_inversion(MGF_ORS_one_hop, p_gamma_th * p_N0/(kappa* p_P0));
    pr_anal_ORS2_outage(Pindex) = euler_inversion(MGF_ORS2_one_hop, p_gamma_th * p_N0/(kappa* p_P0));
    pr_anal_ORS3_outage(Pindex) = euler_inversion(MGF_ORS3_one_hop, p_gamma_th * p_N0/(kappa* p_P0));
    
    % ORSJ
    
    % sum part
    ORSJ_sum = 0;
    ORSJ2_sum = 0;
    ORSJ3_sum = 0;
    
    % q1
    for k = 0:1:1
        for j = 0:1:(1-k)            
            ORSJ_sum =  ORSJ_sum + ((-1)^j * pochhammer(k-1,j))/(factorial(j)) ...
                       * (p_P0/(2*p_Pj*p_gamma_th))^((j+k)/2) ...
                       * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                       * besseli(j+k, 2*V*sqrt(2 * p_P0 * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th));
        end
    end
    
    pr_anal_ORSJ_outage(Pindex) = 1 + 2 * (1-q) * q * ( ...  
                                  exp(-V*(p_P0 + n * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th)) ...  
                                  * p_Pj * p_gamma_th * (1 + V)/(p_P0 + p_Pj * p_gamma_th)...
                                  * besseli(0, 2*V*sqrt(p_P0 * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th))...
                                  - marcumq(sqrt(2*p_P0 * V/(p_P0 + p_Pj * p_gamma_th)), ...
                                            sqrt(2*V*p_Pj * p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 1)) ...
                                  + q^2*(2^((n-1)/2)*exp(-V*(2*p_P0 + n * p_Pj * p_gamma_th)...
                                                         / (p_P0 + p_Pj * p_gamma_th)) ...
                                         * n * p_Pj^2 * p_gamma_th^2/(p_P0 + p_Pj * p_gamma_th)^2 ...
                                         * ORSJ_sum ...
                                         - marcumq(2*sqrt(p_P0 * V/(p_P0 + p_Pj * p_gamma_th)), ...
                                                   sqrt(2 * V * p_Pj * p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), ...
                                                   2))
    % q2                                           
    for k = 0:1:1
        for j = 0:1:(1-k)            
            ORSJ2_sum =  ORSJ2_sum + ((-1)^j * pochhammer(k-1,j))/(factorial(j)) ...
                       * (p_P0/(2*p_Pj*p_gamma_th))^((j+k)/2) ...
                       * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                       * besseli(j+k, 2*V*sqrt(2 * p_P0 * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th));
        end
    end
           
        pr_anal_ORSJ2_outage(Pindex) = 1 + 2 * (1-q2) * q2 * ( ...  
                                  exp(-V*(p_P0 + n * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th)) ...  
                                  * p_Pj * p_gamma_th * (1 + V)/(p_P0 + p_Pj * p_gamma_th)...
                                  * besseli(0, 2*V*sqrt(p_P0 * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th))...
                                  - marcumq(sqrt(2*p_P0 * V/(p_P0 + p_Pj * p_gamma_th)), ...
                                            sqrt(2*V*p_Pj * p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 1)) ...
                                  + q2^2*(2^((n-1)/2)*exp(-V*(2*p_P0 + n * p_Pj * p_gamma_th)...
                                                         / (p_P0 + p_Pj * p_gamma_th)) ...
                                         * n * p_Pj^2 * p_gamma_th^2/(p_P0 + p_Pj * p_gamma_th)^2 ...
                                         * ORSJ_sum ...
                                         - marcumq(2*sqrt(p_P0 * V/(p_P0 + p_Pj * p_gamma_th)), ...
                                                   sqrt(2 * V * p_Pj * p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), ...
                                                   2))
    % q3
    for k = 0:1:1
        for j = 0:1:(1-k)            
            ORSJ3_sum =  ORSJ3_sum + ((-1)^j * pochhammer(k-1,j))/(factorial(j)) ...
                       * (p_P0/(2*p_Pj*p_gamma_th))^((j+k)/2) ...
                       * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                       * besseli(j+k, 2*V*sqrt(2 * p_P0 * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th));
        end
    end
           
        pr_anal_ORSJ3_outage(Pindex) = 1 + 2 * (1-q3) * q3 * ( ...  
                                  exp(-V*(p_P0 + n * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th)) ...  
                                  * p_Pj * p_gamma_th * (1 + V)/(p_P0 + p_Pj * p_gamma_th)...
                                  * besseli(0, 2*V*sqrt(p_P0 * p_Pj * p_gamma_th)/(p_P0 + p_Pj * p_gamma_th))...
                                  - marcumq(sqrt(2*p_P0 * V/(p_P0 + p_Pj * p_gamma_th)), ...
                                            sqrt(2*V*p_Pj * p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 1)) ...
                                  + q3^2*(2^((n-1)/2)*exp(-V*(2*p_P0 + n * p_Pj * p_gamma_th)...
                                                         / (p_P0 + p_Pj * p_gamma_th)) ...
                                         * n * p_Pj^2 * p_gamma_th^2/(p_P0 + p_Pj * p_gamma_th)^2 ...
                                         * ORSJ_sum ...
                                         - marcumq(2*sqrt(p_P0 * V/(p_P0 + p_Pj * p_gamma_th)), ...
                                                   sqrt(2 * V * p_Pj * p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), ...
                                                   2))                                     
                                               
                                               
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
                               
    pr_anal_ORSMJ_outage(Pindex) = 1 ...
                                 - (2 * (1-q) * q * marcumq( sqrt(2*p_P0*V/(p_P0+p_Pj*p_gamma_th)),sqrt(2*n*V*p_Pj*p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 1) ...
                                 - 2* (1-q) * q * (1+V)/sqrt(n * V) * exp(- V - (n-1) * p_Pj * V * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) ...
                                 * (p_P0* (1+V)/(p_Pj * p_gamma_th * Omega* n * V)) ^ ((n-1)/2) ...
                                 * (p_Pj * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) ^ n ...
                                 * (V * Omega / (1+V))^(n/2) ...
                                 * ORSMJ_sum_1 ...
                                 + q^2 * marcumq(2 * sqrt(p_P0 * V / (p_P0 + p_Pj * p_gamma_th)), sqrt(2*n*V*p_Pj*p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 2) ...
                                 - q^2 * p_Pj * p_gamma_th * sqrt(n/V) * exp(-2*V - (n-2)*p_Pj * V * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) / (p_P0 + p_Pj * p_gamma_th) ...
                                 * (p_P0 * (1+V) / (p_Pj * p_gamma_th * Omega * n * V)) ^ ((n-1)/2) ...
                                 * (p_Pj * p_gamma_th / (p_P0 + p_Pj * p_gamma_th))^n ...
                                 * (V * Omega / (1+V))^(n/2) ...
                                 * ORSMJ_sum_2)
                             
    %q2                         
    %sum part 1
    ORSMJ2_sum_1 = 0; %initial
    for k = 0:1:n-1
        for j = 0:1:n-k-1
            ORSMJ2_sum_1 = ORSMJ2_sum_1 ...
                        + (-1)^j*pochhammer(1+k-n,j)/factorial(j)...
                        * (p_P0*(V+1)/(p_Pj*p_gamma_th*Omega))^((j+k)/2) ...
                        * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                        * (n*Omega/(1+V))^((j+k+1)/2) ...
                        * besseli(1+j+k-n,2*V*sqrt(p_Pj * p_P0 * p_gamma_th * n) / (p_P0 + p_Pj * p_gamma_th));
        end
    end
    
    %sum part 2
    ORSMJ2_sum_2 = 0; %initial
    for k = 0:1:n
        for j = 0:1:n-k
            ORSMJ2_sum_2 = ORSMJ2_sum_2 ...
                        + (-1)^j * pochhammer(k-n,j) * 2^((n-j-k-1)/2) / factorial(j) ...
                        * (p_P0 * (1+V)/(p_Pj * p_gamma_th * Omega))^((j+k)/2) ...
                        * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                        * (n * Omega / (V+1))^((j+k-1)/2) ...
                        * besseli(1+j+k-n,2 * V * sqrt(2 * p_Pj * p_P0 * p_gamma_th * n) / (p_P0 + p_Pj * p_gamma_th));
        end
    end
                               
    pr_anal_ORSMJ2_outage(Pindex) = 1 ...
                                 - (2 * (1-q2) * q2 * marcumq( sqrt(2*p_P0*V/(p_P0+p_Pj*p_gamma_th)),sqrt(2*n*V*p_Pj*p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 1) ...
                                 - 2* (1-q2) * q2* (1+V)/sqrt(n * V) * exp(- V - (n-1) * p_Pj * V * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) ...
                                 * (p_P0* (1+V)/(p_Pj * p_gamma_th * Omega* n * V)) ^ ((n-1)/2) ...
                                 * (p_Pj * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) ^ n ...
                                 * (V * Omega / (1+V))^(n/2) ...
                                 * ORSMJ2_sum_1 ...
                                 + q2^2 * marcumq(2 * sqrt(p_P0 * V / (p_P0 + p_Pj * p_gamma_th)), sqrt(2*n*V*p_Pj*p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 2) ...
                                 - q2^2 * p_Pj * p_gamma_th * sqrt(n/V) * exp(-2*V - (n-2)*p_Pj * V * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) / (p_P0 + p_Pj * p_gamma_th) ...
                                 * (p_P0 * (1+V) / (p_Pj * p_gamma_th * Omega * n * V)) ^ ((n-1)/2) ...
                                 * (p_Pj * p_gamma_th / (p_P0 + p_Pj * p_gamma_th))^n ...
                                 * (V * Omega / (1+V))^(n/2) ...
                                 * ORSMJ2_sum_2)
 
    %q3                         
    %sum part 1
    ORSMJ3_sum_1 = 0; %initial
    for k = 0:1:n-1
        for j = 0:1:n-k-1
            ORSMJ3_sum_1 = ORSMJ3_sum_1 ...
                        + (-1)^j*pochhammer(1+k-n,j)/factorial(j)...
                        * (p_P0*(V+1)/(p_Pj*p_gamma_th*Omega))^((j+k)/2) ...
                        * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                        * (n*Omega/(1+V))^((j+k+1)/2) ...
                        * besseli(1+j+k-n,2*V*sqrt(p_Pj * p_P0 * p_gamma_th * n) / (p_P0 + p_Pj * p_gamma_th));
        end
    end
    
    %sum part 2
    ORSMJ3_sum_2 = 0; %initial
    for k = 0:1:n
        for j = 0:1:n-k
            ORSMJ3_sum_2 = ORSMJ3_sum_2 ...
                        + (-1)^j * pochhammer(k-n,j) * 2^((n-j-k-1)/2) / factorial(j) ...
                        * (p_P0 * (1+V)/(p_Pj * p_gamma_th * Omega))^((j+k)/2) ...
                        * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                        * (n * Omega / (V+1))^((j+k-1)/2) ...
                        * besseli(1+j+k-n,2 * V * sqrt(2 * p_Pj * p_P0 * p_gamma_th * n) / (p_P0 + p_Pj * p_gamma_th));
        end
    end
                               
    pr_anal_ORSMJ3_outage(Pindex) = 1 ...
                                 - (2 * (1-q3) * q3 * marcumq( sqrt(2*p_P0*V/(p_P0+p_Pj*p_gamma_th)),sqrt(2*n*V*p_Pj*p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 1) ...
                                 - 2* (1-q3) * q3* (1+V)/sqrt(n * V) * exp(- V - (n-1) * p_Pj * V * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) ...
                                 * (p_P0* (1+V)/(p_Pj * p_gamma_th * Omega* n * V)) ^ ((n-1)/2) ...
                                 * (p_Pj * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) ^ n ...
                                 * (V * Omega / (1+V))^(n/2) ...
                                 * ORSMJ3_sum_1 ...
                                 + q3^2 * marcumq(2 * sqrt(p_P0 * V / (p_P0 + p_Pj * p_gamma_th)), sqrt(2*n*V*p_Pj*p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), 2) ...
                                 - q3^2 * p_Pj * p_gamma_th * sqrt(n/V) * exp(-2*V - (n-2)*p_Pj * V * p_gamma_th / (p_P0 + p_Pj * p_gamma_th)) / (p_P0 + p_Pj * p_gamma_th) ...
                                 * (p_P0 * (1+V) / (p_Pj * p_gamma_th * Omega * n * V)) ^ ((n-1)/2) ...
                                 * (p_Pj * p_gamma_th / (p_P0 + p_Pj * p_gamma_th))^n ...
                                 * (V * Omega / (1+V))^(n/2) ...
                                 * ORSMJ3_sum_2)                                           
                             
    % MRCMJ
    
    %q
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
                     
    pr_anal_MRCMJ_outage(Pindex) = 1 ...
                   - q  ...
                   * marcumq(sqrt(2 * V*p_P0/(p_P0 + p_Pj * p_gamma_th)), sqrt(2 * n * V * p_Pj *p_gamma_th / (p_P0 + p_Pj *p_gamma_th)), 1) ...
                   + q*(1+V)*exp(-(V + (n-1)*p_Pj*V*p_gamma_th/(p_P0 + p_Pj * p_gamma_th))) / (Omega * sqrt(n * V)) ...
                   * (p_P0 * (1+V)/(p_Pj * p_gamma_th * Omega * n * V))^((n-1)/2) ...
                   * (p_Pj * V * p_gamma_th/(p_P0 + p_Pj * p_gamma_th))^(n) ...
                   * (Omega/(V*(1+V)))^(n/2) ...
                   * MRCMJ_sum
               
    %q2
    %sum part
    MRCMJ2_sum = 0;
    for k = 0:1:n-1
        for j = 0:1:n-k-1
            MRCMJ2_sum = MRCMJ2_sum ...
                      + (-1)^j * pochhammer(1+k-n,j)/factorial(j) ...
                      * (p_P0 * (1+V)/(p_Pj * p_gamma_th * Omega))^((j+k)/2) ...
                      * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                      * (n * Omega /(1+V))^((j+k+1)/2) ...
                      * besseli(1+j+k-n, 2*V *sqrt(p_P0*p_Pj*n*p_gamma_th)/(p_P0 + p_Pj * p_gamma_th));
        end
    end
                     
    pr_anal_MRCMJ2_outage(Pindex) = 1 ...
                   - q2  ...
                   * marcumq(sqrt(2 * V*p_P0/(p_P0 + p_Pj * p_gamma_th)), sqrt(2 * n * V * p_Pj *p_gamma_th / (p_P0 + p_Pj *p_gamma_th)), 1) ...
                   + q2*(1+V)*exp(-(V + (n-1)*p_Pj*V*p_gamma_th/(p_P0 + p_Pj * p_gamma_th))) / (Omega * sqrt(n * V)) ...
                   * (p_P0 * (1+V)/(p_Pj * p_gamma_th * Omega * n * V))^((n-1)/2) ...
                   * (p_Pj * V * p_gamma_th/(p_P0 + p_Pj * p_gamma_th))^(n) ...
                   * (Omega/(V*(1+V)))^(n/2) ...
                   * MRCMJ2_sum
               
               
    %q3           
    %sum part
    MRCMJ3_sum = 0;
    for k = 0:1:n-1
        for j = 0:1:n-k-1
            MRCMJ3_sum = MRCMJ3_sum ...
                      + (-1)^j * pochhammer(1+k-n,j)/factorial(j) ...
                      * (p_P0 * (1+V)/(p_Pj * p_gamma_th * Omega))^((j+k)/2) ...
                      * ((p_P0 + p_Pj * p_gamma_th)/p_P0)^k ...
                      * (n * Omega /(1+V))^((j+k+1)/2) ...
                      * besseli(1+j+k-n, 2*V *sqrt(p_P0*p_Pj*n*p_gamma_th)/(p_P0 + p_Pj * p_gamma_th));
        end
    end
                     
    pr_anal_MRCMJ3_outage(Pindex) = 1 ...
                   - q3  ...
                   * marcumq(sqrt(2 * V*p_P0/(p_P0 + p_Pj * p_gamma_th)), sqrt(2 * n * V * p_Pj *p_gamma_th / (p_P0 + p_Pj *p_gamma_th)), 1) ...
                   + q3*(1+V)*exp(-(V + (n-1)*p_Pj*V*p_gamma_th/(p_P0 + p_Pj * p_gamma_th))) / (Omega * sqrt(n * V)) ...
                   * (p_P0 * (1+V)/(p_Pj * p_gamma_th * Omega * n * V))^((n-1)/2) ...
                   * (p_Pj * V * p_gamma_th/(p_P0 + p_Pj * p_gamma_th))^(n) ...
                   * (Omega/(V*(1+V)))^(n/2) ...
                   * MRCMJ3_sum
    
    %% Simulation 
    
    % ORS
    
    
    % Initial parameter
    outage_counter_ORS = 0;
    outage_counter_ORS2 = 0;
    outage_counter_ORS3 = 0;
    outage_counter_ORSJ = 0;
    outage_counter_ORSJ2 = 0;
    outage_counter_ORSJ3 = 0;
    outage_counter_ORSMJ = 0;
    outage_counter_ORSMJ2 = 0;
    outage_counter_ORSMJ3 = 0;
    outage_counter_MRCMJ = 0;
    outage_counter_MRCMJ2 = 0;
    outage_counter_MRCMJ3 = 0;
    
    for N = 1:N_max
        
        % ORS
        
        
        % q
        % Initial parameters
        
        % bernoulli distribution
        B1_ORS = random('Binomial',1,q);
        B2_ORS = random('Binomial',1,q);
        
        
        
        % f1
        f1_k_ORS = (random('rician', ...
                            sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
               
        % f2
        f2_k_ORS = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        gamma_inst_ORS = kappa * p_P0 * (B1_ORS * f1_k_ORS + B2_ORS * f2_k_ORS)/p_N0;
        
        
        if gamma_inst_ORS < p_gamma_th
            outage_counter_ORS = outage_counter_ORS + 1;
        end
        
        % q2
        % Initial parameters
        
        % bernoulli distribution
        B1_ORS2 = random('Binomial',1,q2);
        B2_ORS2 = random('Binomial',1,q2);
        
        
        
        % f1
        f1_k_ORS2 = (random('rician', ...
                            sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
               
        % f2
        f2_k_ORS2 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        gamma_inst_ORS2 = kappa * p_P0 * (B1_ORS2 * f1_k_ORS2 + B2_ORS2*f2_k_ORS2)/p_N0;
        
        
        if gamma_inst_ORS2 < p_gamma_th
            outage_counter_ORS2 = outage_counter_ORS2 + 1;
        end
        
        % q3
        % Initial parameters
        
        % bernoulli distribution
        B1_ORS3 = random('Binomial',1,q3);
        B2_ORS3 = random('Binomial',1,q3);
        
        
        
        % f1
        f1_k_ORS3 = (random('rician', ...
                            sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
               
        % f2
        f2_k_ORS3 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        gamma_inst_ORS3 = kappa * p_P0 * (B1_ORS3 * f1_k_ORS3 + B2_ORS3 * f2_k_ORS3)/p_N0;
        
        
        if gamma_inst_ORS3 < p_gamma_th
            outage_counter_ORS3 = outage_counter_ORS3 + 1;
        end
        
        
        
        % ORSJ
        % q
        % bernoulli distribution
        B1_ORSJ = random('Binomial',1,q);
        B2_ORSJ = random('Binomial',1,q);
        
        
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
        
        % q2
        % bernoulli distribution
        B1_ORSJ2 = random('Binomial',1,q2);
        B2_ORSJ2 = random('Binomial',1,q2);
        
        
        % f1
        f1_ORSJ2 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        % f2
        f2_ORSJ2 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        % v
        v_ORSJ2 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        % gamma_inst
        gamma_inst_ORSJ2 = kappa * p_P0 * (B1_ORSJ2 * f1_ORSJ2 + B2_ORSJ2 * f2_ORSJ2) / (kappa* p_Pj * v_ORSJ2 +p_N0);
        
        if gamma_inst_ORSJ2 < p_gamma_th
            outage_counter_ORSJ2 = outage_counter_ORSJ2 + 1;
        end   
        
        % q3
        % bernoulli distribution
        B1_ORSJ3 = random('Binomial',1,q3);
        B2_ORSJ3 = random('Binomial',1,q3);
        
        
        % f1
        f1_ORSJ3 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        % f2
        f2_ORSJ3 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        % v
        v_ORSJ3 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        % gamma_inst
        gamma_inst_ORSJ3 = kappa * p_P0 * (B1_ORSJ3 * f1_ORSJ3 + B2_ORSJ3 * f2_ORSJ3) / (kappa* p_Pj * v_ORSJ3 +p_N0);
        
        if gamma_inst_ORSJ3 < p_gamma_th
            outage_counter_ORSJ3 = outage_counter_ORSJ3 + 1;
        end   
        
        % ORSMJ
        
         % q
         % Initial parameters
         % bernoulli distribution
        B1_ORSMJ = random('Binomial',1,q);
        B2_ORSMJ = random('Binomial',1,q);
         
        
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
        
         % q2
         % Initial parameters
         % bernoulli distribution
        B1_ORSMJ2 = random('Binomial',1,q2);
        B2_ORSMJ2 = random('Binomial',1,q2);
         
        
        % f1
        f1_ORSMJ2 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
     
        % f2
        f2_ORSMJ2 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        % v_n
        v_n_ORSMJ2 = 0;
        for i = 1:n
            v_n_ORSMJ2 = v_n_ORSMJ2 ...
                      + (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        gamma_inst_ORSMJ2 = kappa * p_P0 * (B1_ORSMJ2 * f1_ORSMJ2 + B2_ORSMJ2 * f2_ORSMJ2) ...
                         / (p_N0 + kappa * p_Pj * v_n_ORSMJ2);
        
        if gamma_inst_ORSMJ2 < p_gamma_th
            outage_counter_ORSMJ2 = outage_counter_ORSMJ2 + 1;
        end
        
         % q3
         % Initial parameters
         % bernoulli distribution
        B1_ORSMJ3 = random('Binomial',1,q3);
        B2_ORSMJ3 = random('Binomial',1,q3);
         
        
        % f1
        f1_ORSMJ3 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
     
        % f2
        f2_ORSMJ3 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        % v_n
        v_n_ORSMJ3 = 0;
        for i = 1:n
            v_n_ORSMJ3 = v_n_ORSMJ3 ...
                       + (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        gamma_inst_ORSMJ3 = kappa * p_P0 * (B1_ORSMJ3 * f1_ORSMJ3 + B2_ORSMJ3 * f2_ORSMJ3) ...
                          / (p_N0 + kappa * p_Pj * v_n_ORSMJ3);
        
        if gamma_inst_ORSMJ3 < p_gamma_th
            outage_counter_ORSMJ3 = outage_counter_ORSMJ3 + 1;
        end
        
        % MRCMJ
        
        
        % q
        % Initial parameters
        % bernoulli distribution
        B1_MRCMJ = random('Binomial',1,q);
        
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
        
        % q2
        % Initial parameters
        % bernoulli distribution
        B1_MRCMJ2 = random('Binomial',1,q2);
        
        % f1
        f1_MRCMJ2 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        % v_n
        v_n_MRCMJ2 = 0;
        for i = 1:n
            v_n_MRCMJ2 = v_n_MRCMJ2 ...
                      + (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        gamma_inst_MRCMJ2 = kappa * p_P0 * B1_MRCMJ2 * f1_MRCMJ2 / (p_N0 + kappa * p_Pj * v_n_MRCMJ2);
        
        if gamma_inst_MRCMJ2 < p_gamma_th
            outage_counter_MRCMJ2 = outage_counter_MRCMJ2 + 1;
        end
        
        % q3
        % Initial parameters
        % bernoulli distribution
        B1_MRCMJ3 = random('Binomial',1,q3);
        
        % f1
        f1_MRCMJ3 = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        % v_n
        v_n_MRCMJ3 = 0;
        for i = 1:n
            v_n_MRCMJ3 = v_n_MRCMJ3 ...
                      + (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        gamma_inst_MRCMJ3 = kappa * p_P0 * B1_MRCMJ3 * f1_MRCMJ3 / (p_N0 + kappa * p_Pj * v_n_MRCMJ3);
        
        if gamma_inst_MRCMJ3 < p_gamma_th
            outage_counter_MRCMJ3 = outage_counter_MRCMJ3 + 1;
        end
        
        
    end
    
    
    pr_simu_ORS_outage(Pindex) = outage_counter_ORS / N_max
    pr_simu_ORS2_outage(Pindex) = outage_counter_ORS2 / N_max
    pr_simu_ORS3_outage(Pindex) = outage_counter_ORS3 / N_max
    
    pr_simu_ORSJ_outage(Pindex) = outage_counter_ORSJ /N_max
    pr_simu_ORSJ2_outage(Pindex) = outage_counter_ORSJ2 /N_max
    pr_simu_ORSJ3_outage(Pindex) = outage_counter_ORSJ3 /N_max
    
    pr_simu_ORSMJ_outage(Pindex) = outage_counter_ORSMJ / N_max
    pr_simu_ORSMJ2_outage(Pindex) = outage_counter_ORSMJ2 / N_max
    pr_simu_ORSMJ3_outage(Pindex) = outage_counter_ORSMJ3 / N_max
    
    pr_simu_MRCMJ_outage(Pindex) = outage_counter_MRCMJ / N_max
    pr_simu_MRCMJ2_outage(Pindex) = outage_counter_MRCMJ2 / N_max
    pr_simu_MRCMJ3_outage(Pindex) = outage_counter_MRCMJ3 / N_max
    
    
end


%% Plot

% ORS
p1=semilogy(P0,pr_anal_ORS_outage,'-');
ax = gca;
ax.FontSize=16;
ax.YLim = [0.05,1];
grid on
p1.Color = 'Red';
p1.LineWidth = 2;
xlabel('Transmission Power (P_0)','FontSize',16);
ylabel('Eavesdropping Channel Transmission Outage Probability','FontSize',14);

hold on;

p2 = semilogy(P0,pr_simu_ORS_outage,'v');
p2.MarkerSize = 10;
p2.Color = 'Red';

p1_2 = semilogy(P0,pr_anal_ORS2_outage,'--');
p1_2.Color = 'Red';
p1_2.LineWidth = 2;

p2_2 = semilogy(P0,pr_simu_ORS2_outage,'v');
p2_2.MarkerSize = 10;
p2_2.Color = 'Red';

p1_3 = semilogy(P0,pr_anal_ORS3_outage,'-.');
p1_3.Color = 'Red';
p1_3.LineWidth = 2;

p2_3 = semilogy(P0,pr_simu_ORS3_outage,'v');
p2_3.MarkerSize = 10;
p2_3.Color = 'Red';

% ORSJ

p3=semilogy(P0,pr_anal_ORSJ_outage,'-');
p3.LineWidth = 2;
p3.Color = 'Blue';

p4 = semilogy(P0,pr_simu_ORSJ_outage,'^');
p4.MarkerSize = 10;
p4.Color = 'Blue';

p3_2=semilogy(P0,pr_anal_ORSJ2_outage,'--');
p3_2.LineWidth = 2;
p3_2.Color = 'Blue';

p4_2 = semilogy(P0,pr_simu_ORSJ2_outage,'^');
p4_2.MarkerSize = 10;
p4_2.Color = 'Blue';

p3_3=semilogy(P0,pr_anal_ORSJ3_outage,'-.');
p3_3.LineWidth = 2;
p3_3.Color = 'Blue';

p4_3 = semilogy(P0,pr_simu_ORSJ3_outage,'^');
p4_3.MarkerSize = 10;
p4_3.Color = 'Blue';


% ORSMJ

p5=semilogy(P0,pr_anal_ORSMJ_outage,'-');
p5.LineWidth = 2;
p5.Color = 'Cyan';

p6 = semilogy(P0,pr_simu_ORSMJ_outage,'o');
p6.MarkerSize = 10;
p6.Color = 'Cyan';

p5_2=semilogy(P0,pr_anal_ORSMJ2_outage,'--');
p5_2.LineWidth = 2;
p5_2.Color = 'Cyan';

p6_2 = semilogy(P0,pr_simu_ORSMJ2_outage,'o');
p6_2.MarkerSize = 10;
p6_2.Color = 'Cyan';

p5_3=semilogy(P0,pr_anal_ORSMJ3_outage,'-.');
p5_3.LineWidth = 2;
p5_3.Color = 'Cyan';

p6_3 = semilogy(P0,pr_simu_ORSMJ3_outage,'o');
p6_3.MarkerSize = 10;
p6_3.Color = 'Cyan';

% MRCMJ

p7=semilogy(P0,pr_anal_MRCMJ_outage,'-');
p7.LineWidth = 2;
p7.Color = 'Magenta';

p8=semilogy(P0,pr_simu_MRCMJ_outage,'+');
p8.MarkerSize = 10;
p8.Color = 'Magenta';

p7_2=semilogy(P0,pr_anal_MRCMJ2_outage,'--');
p7_2.LineWidth = 2;
p7_2.Color = 'Magenta';

p8_2=semilogy(P0,pr_simu_MRCMJ2_outage,'+');
p8_2.MarkerSize = 10;
p8_2.Color = 'Magenta';

p7_3=semilogy(P0,pr_anal_MRCMJ3_outage,'-.');
p7_3.LineWidth = 2;
p7_3.Color = 'Magenta';

p8_3=semilogy(P0,pr_simu_MRCMJ3_outage,'+');
p8_3.MarkerSize = 10;
p8_3.Color = 'Magenta';


% Line

p9=semilogy(P0,pr_line_MRCMJ,'.-.');
p9.LineWidth = 2;
p9.Color = 'Black';

p10=semilogy(P0,pr_line_other,':');
p10.LineWidth = 2;
p10.Color = 'Black';

p9_2=semilogy(P0,pr_line_MRCMJ2,'.-.');
p9_2.LineWidth = 2;
p9_2.Color = 'Black';

p10_2=semilogy(P0,pr_line_other2,':');
p10_2.LineWidth = 2;
p10_2.Color = 'Black';

p9_3=semilogy(P0,pr_line_MRCMJ3,'.-.');
p9_3.LineWidth = 2;
p9_3.Color = 'Black';

p10_3=semilogy(P0,pr_line_other3,':');
p10_3.LineWidth = 2;
p10_3.Color = 'Black';

% Line
p0=semilogy(P0,pr_line_index,'-');
p0.LineWidth = 2;
p0.Color = 'Black';

p0_2 = semilogy(P0,pr_line_index2,'--');
p0_2.LineWidth = 2;
p0_2.Color = 'Black';

p0_3 = semilogy(P0,pr_line_index3,'-.');
p0_3.LineWidth = 2;
p0_3.Color = 'Black';

% Brand
lgd=legend([p0,p0_2,p0_3,p2,p4,p6,p8,p10,p9],'Anal. q=0.3','Anal. q=0.5','Anal. q=0.7','ORS','ORSJ','ORSMJ','MRCMJ','Eq. (66)', 'Eq. (67)');
lgd.Location = 'southwest';
lgd.FontSize = 14;

fname = '/users/shin/dropbox/programming/matlab';
saveas(p1,fullfile(fname,'sub_fig_j_c1'),'fig');
