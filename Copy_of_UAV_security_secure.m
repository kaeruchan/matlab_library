% Clean all residual data
clear all;
close all;
clc;

% parameter library
N_max = 1000;
n_MAX = 20;
M=10; %The number of relay
V=5; %Rician Factor
Omega = 2; 
Pj = 15; %Power of Jammer (dBm)
P0 = 0:20; %Power of Signal (dBm)
gamma_k = 20; %threshold SNR (dBm)
N0 = -150; %Power of noise (dBm)
%beta = 0.5; %time ratio
q = 0.1; % Bernoulli Variance
R = 1; % secrecy rate

% initial parameter
pr_secrecy_analysis = zeros(1,length(P0));
pr_secrecy_simulation = zeros(1,length(P0));

for Pindex = 1:length(P0)
    
     % Parameters (dbm --> watt)    
    p_Pj = 10^(Pj/10)*10^(-3);
    p_P0 = 10^(P0(Pindex)/10)*10^(-3);
    p_N0 = 10^(N0/10)*10^(-3);
    
  

      
      
     pr_secrecy=0;
     f_gamma_e = 0;
    for n=1:n_MAX

      
        % f_gamma_e
    
        sum_T2=@(x) 0;
        sum_T4=@(x) 0;
        %T2
        %sum h,j        
        for k=0:1:M-2
            for j=0:1:M-k-2
                sum_T2 =@(x) sum_T2(x) ...
                       + (-1).^j .* (M-1)*p_Pj ./(p_P0*factorial(j)) ...
                       .* (p_P0*(M-1)./(p_Pj*x)).^((j+k)./2) ...
                       .* ((V+1)*p_P0*p_Pj*x./((M-1).*(p_P0+p_Pj*x).^2)).^(M./2) ...
                       .* (1+p_Pj*x./p_P0).^k ...
                       .* besseli(2+j+k-M,2*V*sqrt((M-1)*p_P0*p_Pj*x)./(p_P0+p_Pj*x)) ...
                       .* (k+1) .* exp(-V*(p_P0+(M-1)*p_Pj*x)./(p_P0+p_Pj*x)) .* pochhammer(2+k-M,j) ... ..//1
                       + (-1).^j .* M .* (M-1).* (p_P0-p_Pj*x)./(2*p_P0*x*factorial(j)) ...
                       .* (p_P0*(M-1)./(p_Pj*x)).^((j+k)./2) .* ((V+1)*p_P0*p_Pj*x./((M-1).*(p_P0+p_Pj*x).^2)).^(M./2) ...
                       .* (1+p_Pj*x./p_P0).^k .* besseli(2+j+k-M,2*V*sqrt((M-1)*p_P0*p_Pj*x)./(p_P0+p_Pj*x)) ...
                       .* exp(-V*(p_P0+(M-1)*p_Pj*x)./(p_P0+p_Pj*x)) ...
                       .* pochhammer(2+k-M,j) ... ..//2
                       + (-1).^(j+1) .* (j+k) .* (1+V) .* p_Pj ./ (2 .* p_P0 .* x .* factorial(j)) ...
                       .* (p_P0*(M-1)./(p_Pj*x)).^((j+k)./2) .* ((V+1)*p_P0*p_Pj*x./((M-1).*(p_P0+p_Pj*x).^2)).^((M-2)./2) ...
                       .* (1+p_Pj*x./p_P0).^k .* besseli(2+j+k-M,2*V*sqrt((M-1)*p_P0*p_Pj*x)./(p_P0+p_Pj*x)) ...
                       .* exp(-V*(p_P0+(M-1)*p_Pj*x)./(p_P0+p_Pj*x)) .* pochhammer(2+k-M,j) ... ..//3
                       + (-1).^(j+1) .* V .* (M-1) .* (M-2) .* p_Pj ./ ( (p_P0 + p_Pj .* x) .* factorial(j) ) ...
                       .* (p_P0*(M-1)./(p_Pj*x)).^((j+k)./2) .* ((V+1)*p_P0*p_Pj*x./((M-1).*(p_P0+p_Pj*x).^2)).^(M./2)  ...
                       .* (1+p_Pj*x./p_P0).^k .* besseli(2+j+k-M,2*V*sqrt((M-1)*p_P0*p_Pj*x)./(p_P0+p_Pj*x)) ...
                       .* exp(-V*(p_P0+(M-1)*p_Pj*x)./(p_P0+p_Pj*x)) .* pochhammer(2+k-M,j) ... ..//4
                       + (-1).^j .* V .* p_Pj .* (M-1).^2 .* (p_P0 - p_Pj .* x) ...
                       ./ (2 .* sqrt((M-1)*p_P0 .* p_Pj *x).*(p_P0 + p_Pj .* x)*factorial(j)) ...
                       .* (p_P0*(M-1)./(p_Pj*x)).^((j+k)./2) .* ((V+1)*p_P0*p_Pj*x./((M-1).*(p_P0+p_Pj*x).^2)).^(M./2)  ...
                       .* (1+p_Pj*x./p_P0).^k .* (besseli(1+j+k-M,2*V*sqrt((M-1)*p_P0*p_Pj*x)./(p_P0+p_Pj*x))  ...
                       + besseli(3+j+k-M,2*V*sqrt((M-1)*p_P0*p_Pj*x)./(p_P0+p_Pj*x)))  ...
                       .* exp(-V*(p_P0+(M-1)*p_Pj*x)./(p_P0+p_Pj*x)) .* pochhammer(2+k-M,j);
            end
        end
        
        T2 =@(x) 2*(1-q)*q ...
           .* (V*(M-1)*p_P0.*p_Pj./(p_P0+p_Pj.*x).^2 ...
           .* exp(-V*(p_P0 + (M-1)*p_Pj*x)./(p_P0+p_Pj.*x)) ...
           .* besseli(0,2.*sqrt(M-1).*V.*sqrt(p_P0.*p_Pj.*x)./(p_P0+p_Pj.*x).^2) ...
           + V*p_P0*p_Pj./(p_P0+p_Pj.*x).^2 ...
           .* (marcumq(sqrt(2*V*p_P0./(p_P0+p_Pj.*x)),sqrt(2*V*(M-1)*p_Pj.*x./(p_P0+p_Pj.*x)),2) ...
           - marcumq(sqrt(2*V*p_P0./(p_P0+p_Pj.*x)),sqrt(2*V*(M-1)*p_Pj.*x./(p_P0+p_Pj*x)),1)) ...
           + sum_T2(x));
        %T4
        %sum k,j
        for k=0:1:M-1
            for j=0:1:M-k-1
               
                sum_T4 =@(x) sum_T4(x) ...
                       + (-1).^j .* 2.^((M-2-j-k)./2) .* (M-1) .* sqrt(V*(M-1))*p_Pj./(p_P0 .* factorial(j)) ...
                       .* (p_P0 .* (M-1)./(p_Pj .* x)).^((j+k)./2) .* (p_P0 .* p_Pj .* x./((M-1).* (p_P0 + p_Pj .* x).^2)).^(M./2) .* (1 + p_Pj .* x./ p_P0).^k .* besseli(2+j+k-M,2 .* V .* sqrt(2 .* (M-1) .* p_P0 .* p_Pj .* x./(p_P0 + p_Pj .* x))) .* exp(-(M + (V+2) .* p_P0) ./ (p_P0 + p_Pj .* x)) .* pochhammer(1+k-M,j) ... ..//1
                       + (-1).^(j+1) .* 2.^((M-j-k-4)./2) .* (j+k) .* (M-1) .* sqrt(V*(M-1)) .* p_Pj ...
                       ./ (p_P0 .* factorial(j)) ...
                       .* (p_P0 .* (M-1)./(p_Pj .* x)).^((j+k)./2) .* (p_P0 .* p_Pj .* x./((M-1).* (p_P0 + p_Pj .* x).^2)).^(M./2) .* (1 + p_Pj .* x./ p_P0).^k .* besseli(2+j+k-M,2 .* V .* sqrt(2 .* (M-1) .* p_P0 .* p_Pj .* x./(p_P0 + p_Pj .* x))) .* exp(-(M + (V+2) .* p_P0) ./ (p_P0 + p_Pj .* x)) .* pochhammer(1+k-M,j) ... ..//2
                       + (-1).^j .* 2.^((M-2-j-k)./2) .* sqrt(V .* (M-1)) .* (M + (V+2) .* p_P0) .* p_Pj.^3 .* x.^2 ...
                       ./ ((p_P0 + p_Pj .* x).^4 .* factorial(j)) ...
                       .* (p_P0 .* (M-1)./(p_Pj .* x)).^((j+k)./2) .* (p_P0 .* p_Pj .* x./((M-1).* (p_P0 + p_Pj .* x).^2)).^((M-2)./2) .* (1 + p_Pj .* x./ p_P0).^k .* besseli(2+j+k-M,2 .* V .* sqrt(2 .* (M-1) .* p_P0 .* p_Pj .* x./(p_P0 + p_Pj .* x))) .* exp(-(M + (V+2) .* p_P0) ./ (p_P0 + p_Pj .* x)) .* pochhammer(1+k-M,j) ... ..//3
                       + (-1).^j .* 2.^((M-j-k-2)./2) .* k .* sqrt(V .* (M-1)) .* p_Pj.^3 .* x.^2 ...
                       ./ ((p_P0 + p_Pj .* x).^3 .* factorial(j)) ... 
                       .* (p_P0 .* (M-1)./(p_Pj .* x)).^((j+k)./2) .* (p_P0 .* p_Pj .* x./((M-1).* (p_P0 + p_Pj .* x).^2)).^((M-2)./2) .* (1 + p_Pj .* x./ p_P0).^k .* besseli(2+j+k-M,2 .* V .* sqrt(2 .* (M-1) .* p_P0 .* p_Pj .* x./(p_P0 + p_Pj .* x))) .* exp(-(M + (V+2) .* p_P0) ./ (p_P0 + p_Pj .* x)) .* pochhammer(1+k-M,j) ... ..//4
                       + (-1).^j .* 2.^((M-j-k-4)./2) .* sqrt(V .* (M-1)) .* M .* p_Pj.^2 .* x ...
                       ./ ((p_P0 + p_Pj .* x).^3 .* factorial(j)) ...
                       .* (p_P0 .* (M-1)./(p_Pj .* x)).^((j+k)./2) .* (p_P0 .* p_Pj .* x./((M-1).* (p_P0 + p_Pj .* x).^2)).^((M-2)./2) .* (1 + p_Pj .* x./ p_P0).^k .* besseli(2+j+k-M,2 .* V .* sqrt(2 .* (M-1) .* p_P0 .* p_Pj .* x./(p_P0 + p_Pj .* x))) .* exp(-(M + (V+2) .* p_P0) ./ (p_P0 + p_Pj .* x)) .* pochhammer(1+k-M,j) ... ..//5
                       + (-1).^j .* 2.^((M-j-k-3)./2) .* exp(-(M + (V+2) .* p_P0) ./ (p_P0 + p_Pj .* x)) .* (M-1).^2 .* (V*(M-1)).^(3/2) .* p_Pj .* pochhammer(1+k-M,j) ...
                       .* sqrt((p_P0 + p_Pj .* x)./((M-1)*p_P0 .* p_Pj .* x)) ...
                       ./ (factorial(j)) ...
                       .*(besseli(1+j+k-M,2 .* V .* sqrt(2 .* (M-1) .* p_P0 .* p_Pj .* x/(p_P0 + p_Pj .* x))) + besseli(3+j+k-M,2 .* V .* sqrt(2 .* (M-1) .* p_P0 .* p_Pj .* x./(p_P0 + p_Pj .* x))));
            end
        end
        
        T4 =@(x) q.^2 .* (V.^2 .* (M-1).^2 .* p_P0 .* p_Pj.^2 .* x ./ (p_P0 + p_Pj .* x).^3 ...
           .* exp(-(V*(2*p_P0 + (M-1)*p_Pj .* x)./(p_P0 + p_Pj .* x))) ...
           .* (p_P0 + p_Pj.* x)./(sqrt(2.*(M-1)).*V.*sqrt(p_P0.*p_Pj.*x)) ...
           .* besseli(1, 2.*sqrt(2).*sqrt(M-1).*V.*sqrt(p_P0.*p_Pj.*x)/(p_P0+p_Pj.*x)) ...
           + 2 .* V .* p_P0 .* p_Pj./(p_P0 + p_Pj .* x).^2 ...
           .* (marcumq(2*sqrt(V*p_P0./(p_P0 + p_Pj .* x)), sqrt(2*V*(M-1)*p_Pj .* x./(p_P0 + p_Pj .* x)),3) ...
           - marcumq(2*sqrt(V*p_P0./(p_P0 + p_Pj .* x)), sqrt(2*V*(M-1)*p_Pj .* x./(p_P0 + p_Pj .* x)),3)) ...
           + sum_T4(x));
    
        f_gamma_e =@(x)  T2(x)+T4(x);


%         f_gamma_e = (2.*exp(-(1+V).*x./Omega + V).*(1-q).*q.*(1+V)./Omega ...
%                   .* besseli(0,2.*sqrt(V.*(1+V).*x./Omega)) ...
%                   + exp(-2 .*V- (V+1).*x ./Omega) .* q^2 ...
%                   .* sqrt((V+1).^3 .*x / (2.*V.*Omega^3)) ...
%                   .* besseli(1,2.*sqrt(2*V.*(1+V).*x./Omega))) ...
%                   .* (1- marcumq(sqrt(2.*(M-1).*V),sqrt(2.*(V+1)/Omega).*sqrt(p_P0.*x - 2.*p_N0)))

    
        % F_gamma_d
    sum_main =@(x) 0;
    %sum part
    for k=0:1:M-2
        for j=0:1:M-k-2
            sum_main =@(x) sum_main(x) + (-1).^(j) .* pochhammer(k+2-M,(j)) ...
                .* (p_P0+p_Pj*x)./(p_P0 .* factorial((j))) ...
                .* exp(-V*(p_P0+(M-1)*p_Pj.*x)./(p_P0+p_Pj.*x)) ...
                .* ((p_P0+p_Pj.*x)./(p_Pj.*x)).^(k) ...
                .* (p_P0./(p_P0+p_Pj.*x)).^M ...
                .* (p_Pj.*x./(p_P0*(M-1))).^(((j)+(k))/2) ...
                .* besseli((j)+(k), 2*V*sqrt(p_P0*(M-1)*p_Pj.*x)./(p_P0+p_Pj.*x));
        end
    end
    
    
    % single channel
        F_gamma_d =@(x) 1- (1-marcumq(sqrt(2*(M-1)*p_Pj*V.*x./(p_P0+p_Pj*x)), ...
                      sqrt(2*V*p_P0./(p_P0+p_Pj.*x)),M-1)-sum_main(x)).^M ...
                      .* (1-marcumq(sqrt(2*(M-1)*p_Pj*V.*x./(p_P0+p_Pj.*x)), ...
                      sqrt(2*V*p_P0./(p_P0+p_Pj.*x)),M-1)-sum_main(x));
    % Totally
       % F_gamma_d =@(x) 1 - (1- F_gamma_d_k(x)).^M..*(1- F_gamma_d_k(x));
        F_gamma =@(x)  F_gamma_d(x) .* f_gamma_e(x);
        pr_secrecy = chebfun(F_gamma,[0,Inf]);      
                    
    end
    pr_secrecy_analysis(Pindex) = pr_secrecy
    
    
    % Simulation
    
    % initial counter
    
    counter = 0;
    
    
    
    for N = 1:N_max
    % gamma_D
        
        
    
    
     % random variables library   
        % jammer
        j_i = 0;
        j_sum_1 = 0;
        j_sum_2 = 0;
        h_2 = 0;
        h_2_max = 0;
        j_2 = 0;
        
        for i = 0:1:M-1
            j_i = j_i + (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
            j_sum_1 = j_sum_1 + (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
            j_sum_2 = j_sum_2 + (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        end
        
         g_2 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
         
        for Mindex = 0:1:M
            h_2 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
            if h_2_max < h_2
                h_2_max = h_2;
            end
        end
        
        for Mindex = 1:M-1
            j_2 = j_2 + (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        end
        
    % Rician fading
        f_1 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;
        f_2 = (random('rician',sqrt(V/(1+V)*Omega),sqrt(Omega/(2*(V+1)))))^2;

    % bernoulli distribution
        i_1 = random('binomial',1,q);
        i_2 = random('binomial',1,q);
        
        
        
        % Gamma_D
        gamma_D = min([(h_2_max*p_P0/(2*p_N0+p_Pj*j_sum_1)), (g_2*p_P0/(2*p_N0+p_Pj*j_sum_2))]);
        gamma_E = (p_P0*i_1*f_1/(2*p_N0+j_2*p_Pj)+p_P0*i_2*f_2/(2*p_N0+j_2*p_Pj));
    
        
        % judge 
        
        C_d = 1/2 .* log2(1+gamma_D);
        
        C_e = 1/2 .* log2(1+gamma_E);
        
        C_s = max(C_d-C_e,0);
        if (C_s < R)
            counter=counter+1;
        end
        
    end
    
    pr_secrecy_simulation(Pindex) = counter/N_max
    
end