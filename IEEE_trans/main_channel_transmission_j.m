%% Clear
clear all;
close all;
clc;

%% Parameter library
N_max = 100000;
M=10; %All relays
K=3; %Number of transmission
n=7; %Number of jammers
Pj = 0:1:20; %Power of Jammer (dBm)
P0 = 30; %Power of Signal (dBm)
gamma_th = 10; %threshold SNR (dB)
N0 = 0; %Power of noise (dBm)
q = 0.5; % Bernoulli Variance
kappa=0.3;

% Channel -- Rician
V = 10; % Rician factor
Omega = 1; % Rician omega


%% Initial Parameters

% ORS
pr_anal_ORS_outage = zeros(1,length(Pj));
pr_simu_ORS_outage = zeros(1,length(Pj));

% ORSJ
pr_anal_ORSJ_outage = zeros(1,length(Pj));
pr_simu_ORSJ_outage = zeros(1,length(Pj));

% ORSMJ

pr_anal_ORSMJ_outage = zeros(1,length(Pj));
pr_simu_ORSMJ_outage = zeros(1,length(Pj));

% MRCMJ

pr_anal_MRCMJ_outage = zeros(1,length(Pj));
pr_simu_MRCMJ_outage = zeros(1,length(Pj));

pr_line_index = zeros(1,length(Pj));
%% Program Start

for Pindex = 1 : length(Pj)
    
    
    %% Parameters  (dBm --> Watt)
    p_Pj = 10^(Pj / 10)*10^(-3); % Jammer
    p_P0 = 10^(P0(Pindex)/10)*10^(-3); % Transmitter
    p_N0 = 10^(N0/10)*10^(-3); % Noise
    p_gamma_th = 10^(gamma_th/10); % threshold SNR
    
    %% Analysis
    
    % ORS
    
    
    % one hop outage
    p_out_one_ORS = 1 - marcumq(sqrt(2*V),sqrt(2*(V+1)/Omega)*sqrt(p_N0*p_gamma_th/(kappa*p_P0)));
    
    pr_anal_ORS_outage(Pindex) = 1 - (1 - (p_out_one_ORS)^M)*(1-p_out_one_ORS);
    
    
    % ORSJ
    
    p_out_one_ORSJ_1 = marcumq(sqrt(2 * V * p_Pj * p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), ...
                       sqrt(2 * V * p_P0 /(p_P0 + p_Pj * p_gamma_th)),1) ...
                     - exp(-V)*p_P0/(p_P0 + p_Pj * p_gamma_th) ...
                     * besseli(0, 2*V*sqrt(p_P0 * p_Pj * p_gamma_th)/(p_P0 + p_Pj*p_gamma_th));
    
    pr_anal_ORSJ_outage(Pindex) = 1 - (1 - (p_out_one_ORSJ)^M)...
                        * (1- p_out_one_ORSJ);
                 
    % ORSMJ
    
    p_out_one_ORSMJ = 1 - marcumq(sqrt(2*V), ...
                                  sqrt(2*(V+1) * p_N0 * p_gamma_th / (kappa * p_P0 * Omega)),1);
    
                              
                              
    pr_anal_ORSMJ_outage(Pindex) = 1 - (1 - (p_out_one_ORSMJ)^M)...
                                 * (1- p_out_one_ORSMJ);
                             
    % MRCMJ
    
    % Rician
%     MGF_MRCMJ_main =@(s) (2 .* (V+1) .* p_N0 * exp(-(V.* p_P0 .* s .* kappa .* Omega)./((V+1)*p_N0 + p_P0 .* s .* kappa .* Omega))...
%                          ./ ((V+1)*p_N0 + p_P0 .* s .* kappa .* Omega)...
%                          + exp(-V .* p_P0 .* s .* kappa .* Omega ...
%                          ./ (2 .* (V + 1) .* p_N0 + p_P0 .* s .* kappa .* Omega)) .* p_P0 .* kappa .* Omega ...
%                          ./ (2 .* (V + 1) .* p_N0 + p_P0 .* s .* kappa .* Omega) ...
%                          .* besseli(0, 2*V*(V+1)*p_N0/(2*(V+1)*p_N0 + p_P0 .* s .* kappa .* Omega)) ...
%                          - p_P0 .* kappa .* Omega .* exp((V+1)*V*p_N0 ./ (p_N0 .* (V+1) + p_P0 .* s .* kappa .* Omega)) ...
%                          ./ (p_N0 .* (V + 1) + p_P0 .* s .* kappa .* Omega) ...
%                          .* MarcumQ(1,...
%                          (V+1)*p_N0 * sqrt(2 * V / ((p_P0 * kappa * Omega * s + (V + 1) * p_N0) * (p_P0 * kappa * Omega * s + 2 * (V + 1) * p_N0))), ...
%                          sqrt(2*V*((V+1)*p_N0 + p_P0 * s * kappa * Omega)...
%                          /(2 * (V+1) * p_N0 + p_P0 * s * kappa * Omega))...
%                          ))^M ./ s;
%                      
%     CDF_MRCMJ_main = euler_inversion(MGF_MRCMJ_main, gamma_th*p_N0/(kappa*p_P0)); 

    % Rician --> Nakagami-m
%     m = (V^2 + 2 * V + 1) / (2*V + 1);
%     
%     MGF_MRCMJ_main =@(s) (2 * gamma(2*m)/(m*gamma(m))*hypergeom([m,2*m],1+m,-1-s*Omega/m))^K/s;
%     
%     pr_anal_MRCMJ_outage(Pindex) = talbot_inversion(MGF_MRCMJ_main, p_gamma_th*p_N0/(kappa*p_P0))/518899.218734917
    
   pr_anal_MRCMJ_outage = [1.00000000000000,0.999999656034523,0.999999581682673,0.999999519412702,0.999999567874363,0.999999557194082,0.999999571476838,0.999999575977228,0.999998416367832,0.999428853442468,0.972726448802497,0.763911617282026,0.355836048958630,0.0876443762896165,0.0122017705329071,0.00107910880678880,6.74989558416799e-05,3.24740417006059e-06,1.27831763414810e-07,4.30761748296529e-09,1.28507335704349e-10,3.48202922409393e-12,8.74174724472822e-14,2.06554799601882e-15,4.65093935657919e-17,1.00786741556083e-18,2.11853111892161e-20,4.34658957011005e-22,8.74788408258316e-24,1.73386060729241e-25,3.39505235056163e-27,6.58391760522107e-29,1.26703848848074e-30,2.42352581000682e-32,4.61319846293843e-34,8.74753489457461e-36];
   
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
        
        % h
        h_ORS_max = 0;
        
        
        for n_count = 1:n
            h_k_ORS = (random('rician', ...
                            sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
            
            if h_k_ORS > h_ORS_max
                h_ORS_max = h_k_ORS;
            end
        end
        
        
        % g
        g_k_ORS = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        
        gamma_inst_s_r_ORS = kappa * p_P0 * h_ORS_max / p_N0;
        
        gamma_inst_r_d_ORS = kappa * p_P0 * g_k_ORS / p_N0;
        
        gamma_inst_ORS = min(gamma_inst_s_r_ORS,gamma_inst_r_d_ORS);
        
        
        if gamma_inst_ORS < p_gamma_th
            outage_counter_ORS = outage_counter_ORS + 1;
        end
        
        
        % ORSJ
        
        % h
        h_ORSJ_max = 0;
        
        
        for n_count_ORSJ = 1:n
            h_k_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
            
            if h_k_ORSJ > h_ORSJ_max
                h_ORSJ_max = h_k_ORSJ;
            end
        end
        
        % g 
        g_k_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        % j
        j_ORSJ_s_r = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        j_ORSJ_r_d = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        gamma_inst_s_r_ORSJ = kappa * p_P0 * h_ORSJ_max / (p_N0 + kappa * p_Pj * j_ORSJ_s_r);
        
        gamma_inst_r_d_ORSJ = kappa * p_P0 * g_k_ORSJ / (p_N0 + kappa * p_Pj * j_ORSJ_r_d);
        
        gamma_inst_ORSJ = min(gamma_inst_s_r_ORSJ,gamma_inst_r_d_ORSJ);
        
        if gamma_inst_ORSJ < p_gamma_th
            outage_counter_ORSJ = outage_counter_ORSJ + 1;
        end   
        
        % ORSMJ
        
         % Initial parameters
        
        % h
        h_ORSMJ_max = 0;
        
        
        for n_count = 1:n
            h_k_ORSMJ = (random('rician', ...
                            sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
            
            if h_k_ORSMJ > h_ORSMJ_max
                h_ORSMJ_max = h_k_ORSMJ;
            end
        end
        
        
        % g
        g_k_ORSMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        
        gamma_inst_s_r_ORSMJ = kappa * p_P0 * h_ORSMJ_max / p_N0;
        
        gamma_inst_r_d_ORSMJ = kappa * p_P0 * g_k_ORSMJ / p_N0;
        
        gamma_inst_ORSMJ = min(gamma_inst_s_r_ORSMJ,gamma_inst_r_d_ORSMJ);
        
        
        if gamma_inst_ORSMJ < p_gamma_th
            outage_counter_ORSMJ = outage_counter_ORSMJ + 1;
        end
        
        % MRCMJ
        
        % one link
        
        % initial parameter
        gamma_MRCMJ = 0;
        
        for kindex = 1 : K
            h_k_MRCMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
            g_k_MRCMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
            
            gamma_s_r_MRCMJ = kappa * p_P0 * h_k_MRCMJ/p_N0;
            gamma_r_d_MRCMJ = kappa * p_P0 * g_k_MRCMJ/p_N0;
            
            gamma_MRCMJ = gamma_MRCMJ ...
                        + min(gamma_s_r_MRCMJ,gamma_r_d_MRCMJ);
        end
        if gamma_MRCMJ < p_gamma_th
            outage_counter_MRCMJ = outage_counter_MRCMJ + 1;
        end
        
        
                
    end
    
    
    pr_simu_ORS_outage(Pindex) = outage_counter_ORS / N_max;
    
    pr_simu_ORSJ_outage(Pindex) = outage_counter_ORSJ /N_max;
    
    pr_simu_ORSMJ_outage(Pindex) = outage_counter_ORSMJ / N_max;
    
    pr_simu_MRCMJ_outage(Pindex) = outage_counter_MRCMJ / N_max
end


%% Plot

% ORS
p1=semilogy(P0,pr_anal_ORS_outage,'-');
ax = gca;
ax.FontSize=16;
ax.YLim = [0.0001,1];
grid on
p1.Color = 'Red';
p1.LineWidth = 2;
xlabel('Transmission Power (P_0)','FontSize',16);
ylabel('Main Channel Transmission Outage Probability','FontSize',16);

hold on;

p2 = semilogy(P0,pr_simu_ORS_outage,'v');
p2.MarkerSize = 10;
p2.Color = 'Red';

% ORSJ

p3=semilogy(P0,pr_anal_ORSJ_outage,'-');
p3.LineWidth = 2;
p3.Color = 'Blue';

p4 = semilogy(P0,pr_simu_ORSJ_outage,'^');
p4.MarkerSize = 10;
p4.Color = 'Blue';

% ORSMJ

p5=semilogy(P0,pr_anal_ORSMJ_outage,'-');
p5.LineWidth = 2;
p5.Color = 'Cyan';

p6=semilogy(P0,pr_simu_ORSMJ_outage,'o');
p6.MarkerSize = 10;
p6.Color = 'Cyan';

% MRCMJ

p7=semilogy(P0,pr_anal_MRCMJ_outage,'-');
p7.LineWidth=2;
p7.Color = 'Magenta';

p8=semilogy(P0,pr_simu_MRCMJ_outage,'+');
p8.MarkerSize=10;
p8.Color = 'Magenta';

% Line
p0=semilogy(P0,pr_line_index,'-');
p0.LineWidth = 2;
p0.Color = 'Black';

% Brand
lgd=legend([p0,p2,p4,p6,p8],'Anal.','ORS','ORSJ','ORSMJ','MRCMJ');
lgd.Location = 'southwest';
lgd.FontSize = 14;

fname = '/users/shin/dropbox/programming/matlab';
saveas(p1,fullfile(fname,'main_fig_j_c2'),'fig');