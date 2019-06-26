%% Clear
clear all;close all;clc;

%% Parameter library
N_max = 50000;
M=10; %All relays
K=5; %Number of transmission
n=3; %Number of jammers
Pj = 15; %Power of Jammer (dBm)
P0 = 0:2:50; %Power of Signal (dBm)
gamma_th = 10; %threshold SNR (dB)
N0 = 0; %Power of noise (dBm)
q = 0.7; % Bernoulli Variance
kappa=0.3;

% Channel -- Rician
V = 10; % Rician factor
Omega = 1; % Rician omega

%% Initial Parameters

pr_anal_ORS_secrecy = zeros(1,length(P0));
pr_simu_ORS_secrecy = zeros(1,length(P0));

pr_anal_ORSJ_secrecy = zeros(1,length(P0));
pr_simu_ORSJ_secrecy = zeros(1,length(P0));

pr_anal_ORSMJ_secrecy = zeros(1,length(P0));
pr_simu_ORSMJ_secrecy = zeros(1,length(P0));

pr_anal_MRCMJ_secrecy = zeros(1,length(P0));
pr_simu_MRCMJ_secrecy = zeros(1,length(P0));

%% Program Part

for Pindex = 1 : length(P0)
    
   % Parameter initial
    p_Pj = 10^(Pj / 10)*10^(-3); % Jammer
    p_P0 = 10^(P0(Pindex)/10)*10^(-3); % Transmitter
    p_N0 = 10^(N0/10)*10^(-3); % Noise
    p_gamma_th = 10^(gamma_th/10); % threshold SNR
    
    %% Analysis
    
    
    % ORS
    
    
    f_e_ORS = @(x) q .* exp(-(x.*(1+V)./Omega +V)) .* ((1+V)/Omega) .* besseli(0,2.*sqrt(V.*(1+V).*x/Omega));
    
    F_d_ORS = @(x) 1 - (1 - (1 - marcumq(sqrt(2.*V),sqrt(2.*(1+V).*(p_gamma_th.*(1+x)-1)/Omega),1)).^K)...
            .* (1 - marcumq(sqrt(2*V),sqrt(2.*(1+V).*(p_gamma_th.*(1+x)-1)/Omega),1));
    
    
    pr_anal_ORS_secrecy(Pindex) = integral(@(x) F_d_ORS(x) .* f_e_ORS(x),0,Inf);
     

    % ORSJ
    
    
    
    %% Simulation 
    
    % initial parameter
    outage_counter_ORS = 0;
    outage_counter_ORSJ = 0;
    outage_counter_ORSMJ = 0;
    outage_counter_MRCMJ = 0;
    
    
    
    
    
    
    
    for N = 1:N_max        
        %% ORS
        
        % main channel
        
        %s_r
        h_ORS_max = 0;
        for kindex=1:K           
            h_ORS = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
            
            if h_ORS > h_ORS_max
                h_ORS_max = h_ORS;
            end
        end
        
        
        g_ORS = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
         
        
        
        gamma_s_r_ORS = kappa * p_P0 * h_ORS_max / p_N0;
        gamma_r_d_ORS = kappa * p_P0 * g_ORS / p_N0;
        
        gamma_main_ORS = min(gamma_s_r_ORS, gamma_r_d_ORS);
        
        
        
        %sub channel
        
        %s_e
        % bernoulli distribution
        B1_ORS = random('Binomial',1,q);
        B2_ORS = random('Binomial',1,q);
        f1_ORS = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        f2_ORS = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        
        gamma_e_ORS = kappa * p_P0 * (B1_ORS * f1_ORS + B2_ORS * f2_ORS) / p_N0;
        
        
        
        if ((1+gamma_main_ORS)/(1+gamma_e_ORS) < p_gamma_th)
           outage_counter_ORS = outage_counter_ORS + 1;
        end
        
        %% ORSJ
        
        % main
        % s-r
        gamma_s_r_ORSJ_max = 0;
        for kindex=1:K
            h_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
            j_1_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
            gamma_s_r_ORSJ = kappa * p_P0 * h_ORSJ / (p_N0 + kappa * p_Pj * j_1_ORSJ);
            if gamma_s_r_ORSJ_max < gamma_s_r_ORSJ
                gamma_s_r_ORSJ_max = gamma_s_r_ORSJ;
            end
        end      
        
        
        % r-d
        g_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        j_2_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        gamma_r_d_ORSJ = kappa * p_P0 * g_ORSJ / (p_N0 + kappa * p_Pj * j_2_ORSJ);
        
        gamma_main_ORSJ = min (gamma_s_r_ORSJ_max,gamma_r_d_ORSJ);
        
        % sub channel
        B1_ORSJ = random('Binomial',1,q);
        B2_ORSJ = random('Binomial',1,q);
        f1_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        f2_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        v_ORSJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        
        
        
        gamma_s_e_ORSJ = kappa * p_P0 * B1_ORSJ * f1_ORSJ / (p_N0 + kappa * p_Pj * v_ORSJ);
        
        gamma_r_e_ORSJ = kappa * p_P0 * B2_ORSJ * f2_ORSJ / (p_N0 + kappa * p_Pj * v_ORSJ);
        
        gamma_e_ORSJ = gamma_s_e_ORSJ + gamma_r_e_ORSJ;
        
        if ((1+gamma_main_ORSJ)/(1+gamma_e_ORSJ) < p_gamma_th)
           outage_counter_ORSJ = outage_counter_ORSJ + 1;
        end
        
        %% ORSMJ
        
        % main
        
        % s-r
        gamma_s_r_ORSMJ_max = 0;
        for kindex=1:K
            h_ORSMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
            
            gamma_s_r_ORSMJ = kappa * p_P0 * h_ORSMJ / p_N0;
            if gamma_s_r_ORSMJ_max < gamma_s_r_ORSMJ
                gamma_s_r_ORSMJ_max = gamma_s_r_ORSMJ;
            end
        end
        
        % r-d
        g_ORSMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        
        gamma_r_d_ORSMJ = kappa * p_P0 * g_ORSMJ / p_N0;
        
        gamma_main_ORSMJ = min(gamma_s_r_ORSMJ_max, gamma_r_d_ORSMJ);
        
        % sub channel
        B1_ORSMJ = random('Binomial',1,q);
        B2_ORSMJ = random('Binomial',1,q);
        f1_ORSMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        f2_ORSMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        v_ORSMJ = 0; % multiple
        for nindex=1:n
            v_ORSMJ =v_ORSMJ + (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        gamma_e_ORSMJ = kappa * p_P0 *(B1_ORSMJ * f1_ORSMJ + B2_ORSMJ * f2_ORSMJ) ...
                      / (p_N0 + kappa);
        
        % secrecy          
        if ((1+gamma_main_ORSMJ)/(1+gamma_e_ORSMJ) < p_gamma_th)
           outage_counter_ORSMJ = outage_counter_ORSMJ + 1;
        end
        
        %% MRCMJ
        
        % main channel
        
        % one link
        
        % initial parameter
        gamma_main_MRCMJ = 0;
        
        for kindex = 1 : K
            h_k_MRCMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
            g_k_MRCMJ = (random('rician', sqrt(V*Omega/(1+V)),sqrt(Omega/(2*(V+1)))))^2;
            
            gamma_s_r_MRCMJ = kappa * p_P0 * h_k_MRCMJ/p_N0;
            gamma_r_d_MRCMJ = kappa * p_P0 * g_k_MRCMJ/p_N0;
            
            gamma_main_MRCMJ = gamma_main_MRCMJ ...
                             + min(gamma_s_r_MRCMJ,gamma_r_d_MRCMJ);
        end
        
        % sub channel
        
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
        
        gamma_e_MRCMJ = kappa * p_P0 * B1_MRCMJ * f1_MRCMJ / (p_N0 + kappa * p_Pj * v_n_MRCMJ);
        
        % secrecy          
        if ((1+gamma_main_MRCMJ)/(1+gamma_e_MRCMJ) < p_gamma_th)
           outage_counter_MRCMJ = outage_counter_MRCMJ + 1;
        end
    end
   
    pr_simu_ORS_secrecy(Pindex) = outage_counter_ORS / N_max;
    
    pr_simu_ORSJ_secrecy(Pindex) = outage_counter_ORSJ / N_max;
    
    pr_simu_ORSMJ_secrecy(Pindex) = outage_counter_ORSMJ / N_max;
    
    pr_simu_MRCMJ_secrecy(Pindex) = outage_counter_MRCMJ / N_max;
end


%% Plot

% ORS
p1=semilogy(P0,pr_anal_ORS_outage,'-');
ax = gca;
ax.FontSize=16;
ax.YLim = [0.01,1];
grid on
p1.Color = 'Red';
p1.LineWidth = 2;
xlabel('Transmission Power (P_0)','FontSize',16);
ylabel('Secrecy Outage Probability','FontSize',16);

hold on;

p2 = semilogy(P0,pr_simu_ORS_outage,'v');
p2.MarkerSize = 10;
p2.Color = 'Red';

% ORSJ

% p3=semilogy(P0,pr_anal_ORSJ_outage,'-');
% p3.LineWidth = 2;
% p3.Color = 'Blue';

p4 = semilogy(P0,pr_simu_ORSJ_outage,'^');
p4.MarkerSize = 10;
p4.Color = 'Blue';

% ORSMJ

% p5=semilogy(P0,pr_anal_ORSMJ_outage,'-');
% p5.LineWidth = 2;
% p5.Color = 'Cyan';

p6 = semilogy(P0,pr_simu_ORSMJ_outage,'o');
p6.MarkerSize = 10;
p6.Color = 'Cyan';

% MRCMJ

% p7=semilogy(P0,pr_anal_MRCMJ_outage,'-');
% p7.LineWidth = 2;
% p7.Color = 'Magenta';

p8=semilogy(P0,pr_simu_MRCMJ_outage,'+');
p8.MarkerSize = 10;
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
saveas(p1,fullfile(fname,'secrecy_fig_j'),'fig');
