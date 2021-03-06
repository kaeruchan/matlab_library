%% Clear
clear all;
close all;
clc;

%% Parameter library
N_max = 1;
M=10; %All relays
K=5; %Number of transmission
n=5; %Number of jammers
Pj = 10; %Power of Jammer (dBm)
P0 = 0:1:35; %Power of Signal (dBm)
gamma_th = 10; %threshold SNR (dB)
N0 = 0; %Power of noise (dBm)
q = 0.5; % Bernoulli Variance
kappa=0.3;

% Channel -- Rician
V = 10; % Rician factor
Omega = 1; % Rician omega


%% Initial Parameters

% ORS
pr_anal_ORS_outage = zeros(1,length(P0));
pr_simu_ORS_outage = zeros(1,length(P0));

% ORSJ
pr_anal_ORSJ_outage = zeros(1,length(P0));
pr_simu_ORSJ_outage = zeros(1,length(P0));

% ORSMJ

pr_anal_ORSMJ_outage = zeros(1,length(P0));
pr_simu_ORSMJ_outage = zeros(1,length(P0));

% MRCMJ

pr_anal_MRCMJ_outage = zeros(1,length(P0));
pr_simu_MRCMJ_outage = zeros(1,length(P0));

pr_line_index = zeros(1,length(P0));
%% Program Start

for Pindex = 1 : length(P0)
    
    
    %% Parameters  (dBm --> Watt)
    p_Pj = 10^(Pj / 10)*10^(-3); % Jammer
    p_P0 = 10^(P0(Pindex)/10)*10^(-3); % Transmitter
    %p_Pj = p_P0; % Jammer power = transmitter power
    p_N0 = 10^(N0/10)*10^(-3); % Noise
    p_gamma_th = 10^(gamma_th/10); % threshold SNR
    
    %% Analysis
    
    % ORS
    
    
    % one hop outage
    p_out_one_ORS = 1 - marcumq(sqrt(2*V),sqrt(2*(V+1)/Omega)*sqrt(p_N0*p_gamma_th/(kappa*p_P0)));
    
    pr_anal_ORS_outage(Pindex) = 1 - (1 - (p_out_one_ORS)^M)*(1-p_out_one_ORS)
    
    
    % ORSJ
    
    p_out_one_ORSJ = marcumq(sqrt(2 * V * p_Pj * p_gamma_th/(p_P0 + p_Pj * p_gamma_th)), ...
                       sqrt(2 * V * p_P0 /(p_P0 + p_Pj * p_gamma_th)),1) ...
                     - exp(-V)*p_P0/(p_P0 + p_Pj * p_gamma_th) ...
                     * besseli(0, 2*V*sqrt(p_P0 * p_Pj * p_gamma_th)/(p_P0 + p_Pj*p_gamma_th));
    
    pr_anal_ORSJ_outage(Pindex) = 1 - (1 - (p_out_one_ORSJ)^M)...
                        * (1- p_out_one_ORSJ)
                 
    % ORSMJ
    
    p_out_one_ORSMJ = 1 - marcumq(sqrt(2*V), ...
                                  sqrt(2*(V+1) * p_N0 * p_gamma_th / (kappa * p_P0 * Omega)),1);
    
                              
                              
    pr_anal_ORSMJ_outage(Pindex) = 1 - (1 - (p_out_one_ORSMJ)^M)...
                                 * (1- p_out_one_ORSMJ)
                             
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
%     pr_anal_MRCMJ_outage(Pindex) = talbot_inversion(MGF_MRCMJ_main, p_gamma_th*p_N0/(kappa*p_P0))/3350719624.03586
%     
    pr_anal_MRCMJ_outage = [1.00000000000000,0.999999886089901,0.999999902308867,0.999999905715436,0.999999919581810,0.999999917575763,0.999999911495254,0.999921757597869,0.982103949886835,0.718422116303974,0.212345702118922,0.0203825425042161,0.000735485105203517,1.23906677675909e-05,1.17514219875651e-07,7.23126583839757e-10,3.20132443385473e-12,1.09901669439904e-14,3.09264246737288e-17,7.44015805106789e-20,1.58096327992858e-22,3.04417085325209e-25,5.42025791554011e-28,9.06890395109136e-31,1.44417500767793e-33,2.21117713191688e-36,3.28148326949621e-39,4.75056270924679e-42,6.74313169262276e-45,9.42277743309487e-48,1.30045278627435e-50,1.77712957457777e-53,2.40954318687367e-56,3.24670157660878e-59,4.35310254749452e-62,5.81362004929301e-65]

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
    
    
    pr_simu_ORS_outage(Pindex) = outage_counter_ORS / N_max
    
    pr_simu_ORSJ_outage(Pindex) = outage_counter_ORSJ /N_max
    
    pr_simu_ORSMJ_outage(Pindex) = outage_counter_ORSMJ / N_max
    
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