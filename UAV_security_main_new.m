clear all;close all;clc;

% parameter library
N_max = 10000;
M=10; %The number of relay
V=10; %Rician Factor
Omega = 1; 
Pj = 5; %Power of Jammer (dBm)
P0 = 15:1:30; %Power of Signal (dBm)
gamma_k = 10; %threshold SNR (dB)
N0 = 0; %Power of noise (dBm)
%beta = 0.5; %time ratio
%q = 0.3; % Bernoulli Variance
lambda = 0.1; %Chernoff-bound parameter


% initial parameters
pr_anal_main = zeros(1,length(P0));
pr_simu_main = zeros(1,length(P0));
pr_simple = zeros(1,length(P0));
pr_anal_main_3 = zeros(1,length(P0));
pr_simu_main_3 = zeros(1,length(P0));
pr_anal_main_6 = zeros(1,length(P0));
pr_simu_main_6 = zeros(1,length(P0));


% pr_anal_main_6 = ... 6 point
%     [0.999999999064279,0.999999826394830,0.999984844612860,0.999410905182280,0.990110731926113,0.928713697467235,0.761342151255679,0.540104694592167,0.344889682873851,0.200358218230779,0.107265821910136,0.0541911167127288,0.0264706561025537,0.0127668919079096,0.00618287344834578,0.00304510583158868];
% pr_simu_main_6 = ... 6 point
%     [1,1,1,0.999900000000000,0.998600000000000,0.977400000000000,0.875600000000000,0.673800000000000,0.443100000000000,0.264900000000000,0.146600000000000,0.0751000000000000,0.0368000000000000,0.0189000000000000,0.00860000000000000,0.00420000000000000];
% pr_anal_main_3 = ... 3 point
%     [0.983132291944771,0.927617387753725,0.792518211027552,0.587151181608773,0.379340104782951,0.221919883743400,0.121142368968366,0.0626891692065762,0.0311994886634327,0.0151911424411997,0.00736105497939199,0.00360199785121840,0.00180011855718509,0.000926165155161085,0.000493101727348844,0.000272431986018007];
% pr_simu_main_3 = ... 3 point 
%     [0.999700000000000,0.992900000000000,0.961400000000000,0.840100000000000,0.631700000000000,0.405300000000000,0.238000000000000,0.125000000000000,0.0649000000000000,0.0344000000000000,0.0165000000000000,0.00870000000000000,0.00320000000000000,0.00200000000000000,0.00110000000000000,0.000300000000000000];
pr_anal_line = zeros(1,length(P0));


for Pindex = 1:length(P0)
    

    %M=3;
    M=3;
    
    
    
    % Parameters (dbm --> watt)    
    p_Pj = 10^(Pj/10)*10^(-3);
    p_P0 = 10^(P0(Pindex)/10)*10^(-3);
    p_gamma_k = 10^(gamma_k/10);
    p_N0 = 10^(N0/10)*10^(-3);
    
    % initial all parameters
    sum = 0;

    for k=0:1:M-2
        for j=0:1:M-k-2
            sum = sum + (-1)^j * pochhammer(k+2-M,j) * (p_P0+p_Pj*p_gamma_k)/(p_P0 * factorial(j)) ...
                * exp(-V*(p_P0+(M-1)*p_Pj*p_gamma_k)/(p_P0+p_Pj*p_gamma_k)) ...
                * ((p_P0+p_Pj*p_gamma_k)/(p_Pj*p_gamma_k))^k ...
                * (p_P0/(p_P0+p_Pj*p_gamma_k))^M ...
                * (p_Pj*p_gamma_k/(p_P0*(M-1)))^((j+k)/2) ...
                * besseli(j+k, 2*V*sqrt(p_P0*(M-1)*p_Pj*p_gamma_k)/(p_P0+p_Pj*p_gamma_k));
        end
    end
    
    % single channel
    pr_anal = marcumq(sqrt(2*(M-1)*p_Pj*V*p_gamma_k/(p_P0+p_Pj*p_gamma_k)), ...
                      sqrt(2*V*p_P0/(p_P0+p_Pj*p_gamma_k)),M-1)-sum;
    
    
    
    pr_anal_main_3(Pindex) = 1-(1-pr_anal^M)*(1-pr_anal)        
     
        %simulation
       outage_counter=0;
    
    for N = 1:N_max
        j_sum_1 = 0;
        j_sum_2 = 0;
        h_2_max = 0;
        
        for Mindex = 0:1:M
            h_2 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
            if h_2_max < h_2
                h_2_max = h_2;
            end
        end
        g_2 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        
        for sum=1:(M-1)
            j_sum_1 = j_sum_1+(random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
            j_sum_2 = j_sum_2+(random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        if min([(h_2_max*p_P0/(2*p_N0+p_Pj*j_sum_1)), (g_2*p_P0/(2*p_N0+p_Pj*j_sum_2))]) < p_gamma_k
            outage_counter =outage_counter+1;
        end
    end
    pr_simu_main_3(Pindex) = outage_counter/N_max
    
    
    
    
        %M=6;
    M=6;
    
    
    
    % Parameters (dbm --> watt)    
    p_Pj = 10^(Pj/10)*10^(-3);
    p_P0 = 10^(P0(Pindex)/10)*10^(-3);
    p_gamma_k = 10^(gamma_k/10);
    p_N0 = 10^(N0/10)*10^(-3);
    
    % initial all parameters
    sum = 0;

    for k=0:1:M-2
        for j=0:1:M-k-2
            sum = sum + (-1)^j * pochhammer(k+2-M,j) * (p_P0+p_Pj*p_gamma_k)/(p_P0 * factorial(j)) ...
                * exp(-V*(p_P0+(M-1)*p_Pj*p_gamma_k)/(p_P0+p_Pj*p_gamma_k)) ...
                * ((p_P0+p_Pj*p_gamma_k)/(p_Pj*p_gamma_k))^k ...
                * (p_P0/(p_P0+p_Pj*p_gamma_k))^M ...
                * (p_Pj*p_gamma_k/(p_P0*(M-1)))^((j+k)/2) ...
                * besseli(j+k, 2*V*sqrt(p_P0*(M-1)*p_Pj*p_gamma_k)/(p_P0+p_Pj*p_gamma_k));
        end
    end
    
    % single channel
    pr_anal = marcumq(sqrt(2*(M-1)*p_Pj*V*p_gamma_k/(p_P0+p_Pj*p_gamma_k)), ...
                      sqrt(2*V*p_P0/(p_P0+p_Pj*p_gamma_k)),M-1)-sum;
    
    
    
    pr_anal_main_6(Pindex) = 1-(1-pr_anal^M)*(1-pr_anal)        
     
        %simulation
       outage_counter=0;
    
    for N = 1:N_max
        j_sum_1 = 0;
        j_sum_2 = 0;
        h_2_max = 0;
        
        for Mindex = 0:1:M
            h_2 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
            if h_2_max < h_2
                h_2_max = h_2;
            end
        end
        g_2 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        
        for sum=1:(M-1)
            j_sum_1 = j_sum_1+(random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
            j_sum_2 = j_sum_2+(random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        if min([(h_2_max*p_P0/(2*p_N0+p_Pj*j_sum_1)), (g_2*p_P0/(2*p_N0+p_Pj*j_sum_2))]) < p_gamma_k
            outage_counter =outage_counter+1;
        end
    end
    pr_simu_main_6(Pindex) = outage_counter/N_max
    
   
    
    
    
    
    
            %M=10;
    M=10;
    
    
    
    % Parameters (dbm --> watt)    
    p_Pj = 10^(Pj/10)*10^(-3);
    p_P0 = 10^(P0(Pindex)/10)*10^(-3);
    p_gamma_k = 10^(gamma_k/10);
    p_N0 = 10^(N0/10)*10^(-3);
    
    % initial all parameters
    sum = 0;

    for k=0:1:M-2
        for j=0:1:M-k-2
            sum = sum + (-1)^j * pochhammer(k+2-M,j) * (p_P0+p_Pj*p_gamma_k)/(p_P0 * factorial(j)) ...
                * exp(-V*(p_P0+(M-1)*p_Pj*p_gamma_k)/(p_P0+p_Pj*p_gamma_k)) ...
                * ((p_P0+p_Pj*p_gamma_k)/(p_Pj*p_gamma_k))^k ...
                * (p_P0/(p_P0+p_Pj*p_gamma_k))^M ...
                * (p_Pj*p_gamma_k/(p_P0*(M-1)))^((j+k)/2) ...
                * besseli(j+k, 2*V*sqrt(p_P0*(M-1)*p_Pj*p_gamma_k)/(p_P0+p_Pj*p_gamma_k));
        end
    end
    
    % single channel
    pr_anal = marcumq(sqrt(2*(M-1)*p_Pj*V*p_gamma_k/(p_P0+p_Pj*p_gamma_k)), ...
                      sqrt(2*V*p_P0/(p_P0+p_Pj*p_gamma_k)),M-1)-sum;
    
    
    
    pr_anal_main(Pindex) = 1-(1-pr_anal^M)*(1-pr_anal)        
     
        %simulation
       outage_counter=0;
    
    for N = 1:N_max
        j_sum_1 = 0;
        j_sum_2 = 0;
        h_2_max = 0;
        
        for Mindex = 0:1:M
            h_2 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
            if h_2_max < h_2
                h_2_max = h_2;
            end
        end
        g_2 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        
        for sum=1:(M-1)
            j_sum_1 = j_sum_1+(random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
            j_sum_2 = j_sum_2+(random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        end
        
        if min([(h_2_max*p_P0/(2*p_N0+p_Pj*j_sum_1)), (g_2*p_P0/(2*p_N0+p_Pj*j_sum_2))]) < p_gamma_k
            outage_counter =outage_counter+1;
        end
    end
    pr_simu_main(Pindex) = outage_counter/N_max
    
    
    
end


%plot
p1=semilogy(P0,pr_anal_main,'-');
ax = gca;
ax.FontSize=16;
ax.YLim = [0.0001,1];
grid on
p1.Color='Blue';
p1.LineWidth=2;
xlabel('Transmission Power (P_0)','FontSize',16)
ylabel('Transmission Outage Probability','FontSize',16)

hold on;
p2=semilogy(P0,pr_simu_main,'^');
p2.Color='Blue';
p2.MarkerSize=10;
p2.LineWidth=2;

p3=semilogy(P0,pr_anal_main_3,'-');
p3.Color='Red';
p3.LineWidth=2;

p4=semilogy(P0,pr_simu_main_3,'v');
p4.Color='Red';
p4.MarkerSize=10;
% 
p5=semilogy(P0,pr_anal_main_6,'-');
p5.Color='Green';
p5.LineWidth=2;

p6=semilogy(P0,pr_simu_main_6,'o');
p6.Color='Green';
p6.MarkerSize=10;

p7=semilogy(P0,pr_anal_line,'-');
p7.Color='Black';
p7.LineWidth=2;

lgd=legend([p7,p4,p6,p2],'Anal.','M = 3','M = 6','M = 10');
lgd.FontSize=16;
lgd.Location='southwest';

% fname = '/users/shin/dropbox/programming/matlab';
% saveas(p1,fullfile(fname,'main_fig'),'fig');