% Clean all residual data
clear all;
close all;
clc;

% parameter library
N_max = 10000;
M=3; %The number of relay
V=10; %Rician Factor
Omega = 1; 
Pj = 5; %Power of Jammer (dBm)
P0 = 25; %Power of Signal (dBm)
gamma_k = 10; %threshold SNR (dB)
N0 = 0; %Power of noise (dBm)
%beta = 0.5; %time ratio
q = 0.1:0.1:0.9; % Bernoulli Variance
lambda = 0.1; %Chernoff-bound parameter

% Analysis

pr_anal_sub = zeros(1,length(q));
pr_simu_sub = zeros(1,length(q));
pr_simple = zeros(1,length(q));
pr_anal_sub_infty_3 = zeros(1,length(q));
pr_anal_sub_infty_5 = zeros(1,length(q));
pr_anal_sub_infty_7 = zeros(1,length(q));


%Library (Pr_type_channel_q_N0)
pr_anal_sub_10 = [0.890182450901103,0.782976867476473,0.678383249726110,0.576401597650013,0.477031911248182,0.380274190520618,0.286128435467321,0.194594646088290,0.105672822383526];

pr_simu_sub_10 = [0.907100000000000,0.802200000000000,0.695200000000000,0.606900000000000,0.511800000000000,0.412600000000000,0.300900000000000,0.210100000000000,0.122200000000000];

pr_anal_sub_6 = [0.829314801982305,0.674353423063805,0.535115863244500,0.411602122524389,0.303812200903473,0.211746098381751,0.135403814959224,0.0747853506358914,0.0298907054117535];

pr_simu_sub_6 = [0.840900000000000,0.680200000000000,0.560800000000000,0.422400000000000,0.333800000000000,0.231700000000000,0.150200000000000,0.0903000000000000,0.0387000000000000];

for qindex = 1:length(q)
    
     % Parameters (dbm --> watt)    
    p_Pj = 10^(Pj/10)*10^(-3);
    p_P0 = 10^(P0/10)*10^(-3);
    p_gamma_k = 10^(gamma_k/10);
    p_N0 = 10^(N0/10)*10^(-3);
    
    
    % main
    tau1 = 2*(1-q(qindex))*q(qindex)
    
   % tau3

    tau3 = q(qindex)^2
    

    
    %tau2

    %initial sum part
    sum_tau2=0;
    %sum part
    for k = 0:1:M-2
        for j = 0:1:M-k-2
            
            sum_tau2 = sum_tau2+ (-1)^j *pochhammer(2+k-M,j) /factorial(j) ...
                     * (p_P0*(M-1)/(p_Pj*p_gamma_k))^((j+k)/2) ...
                     * ((p_P0+p_Pj*p_gamma_k)/p_P0)^k ...
                     * besseli(2+j+k-M, 2*V*sqrt(p_P0*p_Pj*(M-1)*p_gamma_k)/(p_P0+p_Pj*p_gamma_k));
            
        end
    end
                 
    tau2 = 2*(1-q(qindex))*q(qindex)*(1 ...
         + exp(-V*(p_P0+(M-1)*p_Pj*p_gamma_k)/(p_P0+p_Pj*p_gamma_k))*(M-1)*(p_P0+p_Pj*p_gamma_k)/p_P0 ...
         * (p_P0*p_Pj*p_gamma_k/((M-1)*(p_P0+p_Pj*p_gamma_k)^2))^(M/2) ...
         * sum_tau2 ...
         - marcumq(sqrt(2*p_P0*V/(p_P0+p_Pj*p_gamma_k)), ...
           sqrt(2*V*(M-1)*p_Pj*p_gamma_k/(p_P0+p_Pj*p_gamma_k)),1))
    
    % tau4
    

    % initial sum part
    sum_tau4=0;
    %sum part
    for k = 0:1:M-1
        for j = 0:1:M-k-1
                  
            sum_tau4 = sum_tau4 ...
                     + (-1)^j * pochhammer(1+k-M,j) / factorial(j) ...
                     * (p_P0*(M-1)/(2*p_Pj*p_gamma_k))^((j+k)/2) ...
                     * ((p_P0 +p_Pj*p_gamma_k)/p_P0)^k ...
                     * besseli(2+j+k-M, ...
                            2*V*sqrt(2*p_P0*p_Pj*(M-1)*p_gamma_k)/(p_P0+p_Pj*p_gamma_k));
        end
    end
    
    tau4 = q(qindex)^2 ...
           * (1 ...
           + 2^((M-2)/2) ...
           * exp(-(V*(2*p_P0+(M-1)*p_Pj*p_gamma_k)/(p_P0+p_Pj*p_gamma_k))) ...
           * (M-1)*p_Pj*p_gamma_k/p_P0 ...
           * (p_P0*p_Pj*p_gamma_k/((M-1)*(p_P0+p_Pj*p_gamma_k)^2))^(M/2) ...
           * sum_tau4 ...
           - marcumq(2*sqrt(p_P0*V/(p_P0+p_Pj*p_gamma_k)), ...
           sqrt(2*V*(M-1)*p_Pj*p_gamma_k/(p_P0+p_Pj*p_gamma_k)),2))
    
    %analysis result
    pr_anal_sub(qindex) = 1-(tau1-tau2+tau3-tau4)
    
        %simulation part
    
    %initial parameters
    outage_counter=0;
    
    for N = 1:N_max
        i_sum_1 = 0;
        i_sum_2 = 0;
        j_2 = 0;
        
        
        for Mindex = 1:M-1
            j_2 = j_2 + (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        
        
        f_1 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        f_2 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;

       %bernoulli distribution
        i_1 = random('binomial',1,q(qindex));
        i_2 = random('binomial',1,q(qindex));
       
        
        if (p_P0*i_1*f_1/(2*p_N0+j_2*p_Pj)+p_P0*i_2*f_2/(2*p_N0+j_2*p_Pj))< p_gamma_k
            outage_counter =outage_counter+1;
        end
    end
    pr_simu_sub(qindex) = outage_counter/N_max
                     
%     pr_anal_sub_infty_3(Pindex) = (1-0.3)^2;   
%     pr_anal_sub_infty_5(Pindex) = (1-0.5)^2;
%     pr_anal_sub_infty_7(Pindex) = (1-0.7)^2;
    
    
end




%plot 
p1=semilogy(q,pr_anal_sub,'-');
ax = gca;
ax.FontSize=16;
ax.YLim = [0.01,1];
grid on
p1.Color='Red';
p1.LineWidth=2;
xlabel('q','FontSize',16)
ylabel('Eavesdropping Channel Transmission Outage Probability','FontSize',12)

hold on;
p2=semilogy(q,pr_simu_sub,'v');
p2.Color='Red';
p2.MarkerSize=10;
p2.LineWidth=2;

p3=semilogy(q,pr_anal_sub_6,'-');
p3.Color='Green';
p3.LineWidth=2;

p4=semilogy(q,pr_simu_sub_6,'o');
p4.Color='Green';
p4.MarkerSize=10;



p17=semilogy(q,pr_anal_sub_10,'-');
p17.Color='Blue';
p17.LineWidth=2;

p18=semilogy(q,pr_simu_sub_10,'^');
p18.Color='Blue';
p18.MarkerSize=10;


p0=semilogy(q,pr_simple,'-');
p0.Color='Black';
p0.LineWidth=2;

lgd=legend([p0,p2,p4,p6],'Anal.','M=3','M=6','M=10');
% set(lgd,'Interpreter','latex')
lgd.FontSize=14;
lgd.Location='southwest';
