clear, close all, clc;


% parameter library
%N=10; %The number of relay
V=5; %Rician Factor
Omega = 1; 
Pj = 20; %Power of Jammer (dBm)
P0 = 0:1:30; %Power of Signal (dBm)
gamma_th = 0; %threshold SNR (dBm)
N0 = -200; %Power of noise (dBm)
% beta = 0.5; %time ratio
q = 0.8; %Bernoulli Variance
M = 3; % Relay number
R =1; %secrecy outage rate
N =10000;
%x_k = 


% Gauss-Laguerre Approximation
% 10point approximation
x_i = [0.137793470540492431,...
    0.729454549503170498,...
    0.180834290174031605e01,...
    0.340143369785489951e01,...
    0.555249614006380363e01,...
    0.833015274676449670e01,...
    0.118437858379000656e02,...
    0.162792578313781021e02,...
    0.219965858119807620e02,...
    0.299206970122738916e02];


w_i = [0.308441115765020141,...
    0.401119929155273552,...
    0.218068287611809422,...
    0.620874560986777475e-01,...
    0.950151697518110055e-02,...
    0.753008388587538775e-03,...
    0.282592334959956557e-04,...
    0.424931398496268637e-06,...
    0.183956482397963078e-08,...
    0.991182721960900856e-12];


% Parameter Initial
anal_outage_m = zeros(1,length(P0));
anal_outage_s = zeros(1,length(P0));
simu_outage_m = zeros(1,length(P0));
simu_outage_s = zeros(1,length(P0));


% analysis

for pIndex = 1:length(P0)
    
    
    %Pj = P0(pIndex);
    % initial for all parameters

    
    
    % Theory part
    %Power transform (dB -> Watt)
    p_Pj = 10^(Pj/10)*10^(-3);
    p_P0 = 10^(P0(pIndex)/10)*10^(-3);
    %p_P0 =;
    p_gamma_th = 10^(gamma_th/10)*10^(-3);
    p_N0 = 10^(N0/10)*10^(-3);

    %sum of main channel initial
    sum_Prm=0;
    B1 = sqrt(2*V);
    B2 = sqrt(2*(1+V)/Omega)*sqrt((2*p_N0+p_Pj)*p_gamma_th/p_P0);
    B3 = 2 * sqrt((V+1)*V*(M-1)/Omega);
    B4 = (V+1)/Omega;
    for k=0:M-2
        for j=0:M-k-2
            sum_Prm = (B1/B3)^(j+k)*((2*B4+B2^2)/B2)^k*(-B2)^j ...
                      *pochhammer(-M+k+2,j)/factorial(j)*besseli(k+j,B3*B1*B2/(2*B4+B2^2)) ...                      
                      +sum_Prm
        end
    end
    sum_Prm = (B3/2)^(M-2)*B4^(-(M-1))*exp(B3^2/(4*B4)) ...
              +2/B3*(B3/(2*B4+B2^2))^(M-1)*exp((B3^2-2*B4*B1)/(2*(2*B4+B2^2))) ...
              *sum_Prm ...
              -1/B4 *(B3/(2*B4))^(M-2)*exp(B3^2/(4*B4)) ...
              *marcumq(B2*B3/sqrt(2*B4*(2*B4+B2^2)),B1*sqrt(2*B4/(2*B4+B2^2)),M-1)
    
    %sum of sub channel initial
    sum_Prs =0;

    % main channel information outage
    Pr_m =  ((V+1)/Omega)^(M/2)*(V*(M-1))^((2-M)/2)*exp(-(M-1)*V) ...
           *sum_Prm
    Pr_m_a = 1 -(1-Pr_m)^M*(1-Pr_m)

    % sub channel information outage
    %Pr_s = 2*q^2*Omega^2 + 2*(1-q)*q*Omega - sum_Prs;
    Pr_s = sum_Prs;
    Pr_s_a = 1 - Pr_s;

    
    anal_outage_m(pIndex) = Pr_m_a;
    anal_outage_s(pIndex) = Pr_s_a;

    % secrecy outage probability
    % Gauss-Laguerre Formula
    
    P_s_out =0; %Initial parameter
   







    % Simulation
    
    %Parameters for simulation
        outage_counter_m =0;
        outage_counter_s =0;
    %Main channel outage
 
    for N_index=1:N
        h_m = zeros(1,M); %Initial parameter

        % Optimal S--R parameter
        h_m_max = random('Rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(1+V))));
% %         for k=1:M
% %             h_m(k) = random('Rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(1+ ...
% %                                                               V))));
% %             if h_m_max < h_m(k)
% %                 h_m_max = h_m(k);
% %             end 
% %         end
        % R--D parameter
        g_m = random('Rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(1+V))));
        
        % S--E parameter
        f_1 = random('Rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(1+V))));
        
        % R*--E parameter
        f_2 = random('Rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(1+V))));
        
        
        % Jammers
        sum_j_i_1 =0; %Initial parameter
        sum_j_i_2 =0;
        sum_v_i_2 =0;
        for i=1:(M-1)
            %Jammer--R
            sum_j_i_1 = sum_j_i_1 + (random('Rician', ...
                                          sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(1+V)))))^2;
            %Jammer--D
            sum_j_i_2 = sum_j_i_2 + (random('Rician', ...
                                          sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(1+V)))))^2;
            %Jammer--E
            sum_v_i_2 = sum_v_i_2 + (random('Rician', ...
                                          sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(1+V)))))^2;
        end
        
        % SINRs --main channel
        SNR_s_r = p_P0*h_m_max^2/(2*p_N0+p_Pj*sum_j_i_1);
        SNR_r_d = p_P0*g_m^2/(2*p_N0+p_Pj*sum_j_i_2);
        
       SNR_m = SNR_s_r;
        
        if SNR_m < p_gamma_th
            outage_counter_m = outage_counter_m + 1;
        end
        
        B_1 = random('Binomial',1,q);
        B_2 = random('Binomial',1,q);
        
        % SNRs -- sub channel
        SNR_s_e = p_P0 * B_1 * (f_1)^2/(2*p_N0 + p_Pj*sum_v_i_2);
        SNR_r_e = p_P0 * B_2 * (f_2)^2/(2*p_N0 + p_Pj*sum_v_i_2);
        
        SNR_b = SNR_s_e + SNR_r_e;
        
        if SNR_b < p_gamma_th
            outage_counter_s =outage_counter_s + 1;
        end

    
    
    
    
    end
    
    simu_outage_m(pIndex) = outage_counter_m/N
    
    simu_outage_s(pIndex) = outage_counter_s/N
end



%Plotting
clear title xlabel ylable
p1=semilogy(P0,simu_outage_m,'-');
ax =gca;
ax.FontSize =16;
ax.YLim = [0.0001,1];
p1.Color='Red';
p1.LineWidth = 2;
grid on
xlabel('P_0 (dBm)','FontSize',16)
ylabel('Outage Probability','FontSize',16)


hold on

p2=semilogy(P0,simu_outage_s,'-');
p2.Color='Blue';
p2.LineWidth = 2;

hold off

lgd = legend([p1,p2],'Main Channel Outage (Simu.)','Sub Channel Outage (Simu.)');
lgd.FontSize =16;
lgd.Location = 'NorthEast';