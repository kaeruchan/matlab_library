clear all;close all;clc;

% parameter library
N_max = 10000;
M=10; %The number of relay
V=5; %Rician Factor
Omega = 1; 
Pj = 5; %Power of Jammer (dBm)
P0 = 0:20; %Power of Signal (dBm)
gamma_k = 20; %threshold SNR (dBm)
N0 = 0; %Power of noise (dBm)
%beta = 0.5; %time ratio
q = 0.3; % Bernoulli Variance
lambda = 0.1; %Chernoff-bound parameter


% initial parameters
pr_anal_main = zeros(1,length(P0));
pr_simu_main = zeros(1,length(P0));
pr_simple = zeros(1,length(P0));

 pr_anal_main_6 = ... 6 point
    [0.869215012140049,0.720286866775624,0.553894072703202,0.404389804613247,0.282940506936184,0.190984782452005,0.125651455757490,0.0814959037078576,0.0526365880465861,0.0341290438186589,0.0223482099803852,0.0148397113504132,0.0100175012693968,0.00688279313916917,0.00481411935607545,0.00342577623504836,0.00247743481892482,0.00181807679430113,0.00135170167954723,0.00101644166478998,0.000771817633713590];
pr_simu_main_6 = ... 6 point
    [0.932400000000000,0.820000000000000,0.649700000000000,0.487100000000000,0.339600000000000,0.229000000000000,0.154600000000000,0.0998000000000000,0.0651000000000000,0.0433000000000000,0.0260000000000000,0.0200000000000000,0.0135000000000000,0.00740000000000000,0.00640000000000000,0.00390000000000000,0.00280000000000000,0.00210000000000000,0.00270000000000000,0.00160000000000000,0.00160000000000000];
pr_anal_main_3 = ... 3 point
    [0.305954505145345,0.206199724957163,0.135976254301641,0.0885113520564509,0.0572896862202077,0.0371237562102699,0.0242284939648363,0.0160012070267748,0.0107289468171602,0.00731757061450422,0.00508039123930415,0.00358970198432163,0.00257908924330430,0.00188166342945240,0.00139186182688245,0.00104207205276285,0.000788358267362899,0.000601717339210550,0.000462685112403771,0.000357977473285254,0.000278374936501735];
pr_simu_main_3 = ... 3 point 
   [0.467600000000000,0.329000000000000,0.214500000000000,0.143500000000000,0.0956000000000000,0.0546000000000000,0.0356000000000000,0.0253000000000000,0.0164000000000000,0.00960000000000000,0.00780000000000000,0.00580000000000000,0.00400000000000000,0.00300000000000000,0.00220000000000000,0.00210000000000000,0.000700000000000000,0.00140000000000000,0.000800000000000000,0.000400000000000000,0.000300000000000000];
pr_anal_line = zeros(1,length(P0));


for Pindex = 1:length(P0)
    

    % Parameters (dbm --> watt)    
    p_Pj = 10^(Pj/10)*10^(-3);
    p_P0 = 10^(P0(Pindex)/10)*10^(-3);
    p_gamma_k = 10^(gamma_k/10)*10^(-3);
    p_N0 = 10^(N0/10)*10^(-3);
    
    % initial all parameters
    sum = 0;
%     alpha = sqrt(2*V);
%     beta = sqrt(2*(1+V)*p_Pj*p_gamma_k/(p_P0*Omega));
%     c = 2*sqrt((1+V)*V*(M-1)/Omega);
%     p = (1+V)/Omega;
%     tilde_p = 2*p+beta^2;
%     mu = M-1;
    
    %sum part
%     for k=0:1:mu-1
%         for j=0:1:mu-k-1
%             sum = sum + (alpha/c)^(j+k) * (tilde_p/beta)^k * (-beta)^j ...
%                 * pochhammer(-mu+k+1 , j)/factorial(j)...
%                 * besseli(k+j,c*alpha*beta/tilde_p);
%         end
%     end

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
%         pr_anal = 1 ...
%                          - ((V+1)/Omega)^(M/2) * (V*(M-1))^(-(M-2)/2)*exp(-(M-1)*V) ...
%                          * ((c/2)^(mu-1)*p^(-mu)*exp(c^2/(4*p))...
%                          + 2/c*(c/tilde_p)^(mu)*exp((c^2/2-p*alpha^2)/tilde_p)...
%                          * sum ...
%                          - 1/p * (c/(2*p))^(mu-1)*exp(c^2/(4*p))...
%                          * marcumq(beta*c/sqrt(2*p*tilde_p),alpha*sqrt(2*p/tilde_p),mu));
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

fname = '/users/shin/dropbox/programming/matlab';
saveas(p1,fullfile(fname,'main_fig'),'fig');