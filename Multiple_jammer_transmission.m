clear all;close all;clc;

% parameter library
N=10; %The number of relay
V=5; %Rician Factor
Omega = 1; 
Pj = 10; %Power of Jammer (dBm)
%P0 = 20; %Power of Signal (dBm)
gamma_k = 20; %threshold SNR (dBm)
N0 = -150; %Power of noise (dBm)
%beta = 0.5; %time ratio
q = 0.2; %Bernoulli Variance
x_k = 




for P0 = 1:20

    
    
    
    
    % initial for all parameters

    
    %infinite sum of the sub-channel 
    L_infinity=10;
    
    % Theory part
    %Power transform (dBm -> Watt)
    p_Pj = 10^(Pj/10)*10^(-3);
    p_P0 = 10^(P0/10)*10^(-3);
    p_gamma_k = 10^(gamma_k/10)*10^(-3);
    p_N0 = 10^(N0/10)*10^(-3);

    %sum of main channel initial
    sum_Prm=0;
    for k = 0:N-1
        for j = 0:N-k-1
            %%%%
        sum_Prm = (-1)^j*pochhammer(k-N+1,j)/factorial(j)*2^(-(2*k+j)/2)...
                  *(p_gamma_k*(2*p_N0+p_Pj)/(Omega*p_P0))^((j-k)/2)...
                  *(V*(1+V)/(N*Omega))^((j+k)/2)...
                  *((1+V)*(p_gamma_k*(2*p_N0+p_Pj))/(p_P0*Omega))^(k/2)...
                  *besseli(j+k,2*sqrt(2)*V*sqrt(N*(1+V)...
                  *(2*p_N0+p_Pj)*p_gamma_k/(Omega*(p_gamma_k*(2*p_N0+p_Pj)-p_P0))))...
                  +sum_Prm
        end
    end

    %sum of sub channel initial
    sum_Prs=0;
    for l = 0 : L_infinity
        for k = 0 : (N+l-1)
            for n = 0 : k
                sum_Prs = sum_Prs+...
                    nchoosek(k,n)...
                    *1/(sqrt(V)*(p_P0+p_Pj*p_gamma_k)*factorial(k)*factorial(l))...
                    *(-1)^n*2^(n-2)...
                    *exp(-(N*V+(4*p_Pj*V*p_gamma_k)/(p_P0+p_Pj*p_gamma_k)-(2*p_N0*(1+V))/(p_Pj*Omega)))...
                    *(N*V)^l*(p_N0/p_Pj)^n*(p_P0/(p_Pj*p_gamma_k))^(k-n)*(1+V/Omega)^k...
                    *((1+V)*(p_P0+p_Pj*p_gamma_k)/(p_Pj*p_gamma_k*Omega))^(n-k-1/2)...
                    *sqrt(Omega/(V*(1+V)))*(2*exp(2*p_Pj*V*p_gamma_k/(p_P0+p_Pj*p_gamma_k))*sqrt(V)...
                    *(p_P0+p_Pj*p_gamma_k)*gamma(1+k-n)...
                    *whittakerM(k-n+1/2,0,4*p_Pj*V*p_gamma_k/(p_P0+p_Pj*p_gamma_k))...
                    +exp(V)*p_Pj*q^2*sqrt(1+V)*p_gamma_k*sqrt((1+V)*(p_P0+p_Pj*p_gamma_k)/(p_Pj*p_gamma_k*Omega))...
                    *Omega^(3/2)*gamma(k-n+2)...
                    *whittakerM(k-n+1,1/2,8*p_Pj*V*p_gamma_k/(p_P0+p_Pj*p_gamma_k)))                          
            end
        end
    end

    %main channel information outage
    Prm = (-1)^N*exp(-N*V)*(N*V)^((N-1)/2)*((1+V)/2)^(-(N+1)/2)...
        *(1-marcumq(...
        sqrt(2*N*V)*sqrt((2*p_N0+p_Pj)*p_gamma_k/(p_P0-(2*p_N0+p_Pj)*p_gamma_k)),...
        sqrt(2*V)*sqrt(p_P0/(p_P0-(2*p_N0+p_Pj)*p_gamma_k)),N))...
        +2^(N/2)*exp(N*(N+2*V)*sqrt((1+V)/(2*(2*N*V-1)*Omega)))*(N*V*(1+V)/Omega)^((N-1)/2)...
        *(((1+V)*(p_gamma_k*(2*p_N0+p_Pj)-p_P0))/(p_P0*Omega))^(-N)...
        +sum_Prm

    Pr_m_a = 1 - Prm*Prm;

    %sub channel information outage
    sum_l=0;

    Pr_s_a = (q^2 * Omega^2 * exp(-V/2)+4*(1-q)*q*Omega)/2 -sum_l


    % Simulation 
    
    
    
end