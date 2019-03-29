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
P0 = 0:50; %Power of Signal (dBm)
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
    p_gamma_k = 10^(gamma_k/10)*10^(-3);
    p_N0 = 10^(N0/10)*10^(-3);
    
    

    % Gauss-Laguerre library (20-point)


    x_k = [0.705398896919887533667e-01, ... 1
           0.372126818001611443794,     ... 2
           0.916582102483273564668, 	... 3
           0.170730653102834388069e01,  ... 4
           0.274919925530943212965e01,  ... 5
           0.404892531385088692237e01,  ... 6
           0.561517497086161651410e01,  ... 7
           0.745901745367106330977e01,  ... 8
           0.959439286958109677247e01,  ... 9
           0.120388025469643163096e02,  ... 10
           0.148142934426307399785e02,  ... 11
           0.179488955205193760174e02,  ... 12
           0.214787882402850109757e02,  ... 13
           0.254517027931869055035e02,  ... 14
           0.299325546317006120067e02,  ... 15
           0.350134342404790000063e02,  ... 16
           0.408330570567285710620e02,  ... 17
           0.476199940473465021399e02,  ... 18
           0.558107957500638988908e02,  ... 19
           0.665244165256157538186e02   ... 20
          ];

    w_k = [0.168746801851113862149,     ... 1
           0.291254362006068281717,     ... 2
           0.266686102867001288550,     ... 3
           0.166002453269506840031,	    ... 4
           0.748260646687923705401e-01, ... 5
           0.249644173092832210728e-01, ... 6
           0.620255084457223684745e-02, ... 7
           0.114496238647690824204e-02, ... 8
           0.155741773027811974780e-03, ... 9
           0.154014408652249156894e-04, ... 10
           0.108648636651798235148e-05, ... 11
           0.533012090955671475093e-07, ... 12
           0.175798117905058200358e-08, ... 13
           0.372550240251232087263e-10, ... 14
           0.476752925157819052449e-12, ... 15
           0.337284424336243841237e-14, ... 16
           0.115501433950039883096e-16, ... 17
           0.153952214058234355346e-19, ... 18
           0.528644272556915782880e-23, ... 19
           0.165645661249902329591e-27  ... 20
          ];


      
      
     pr_secrecy=0;
     f_gamma_e = 0;
     
    for n=1:n_MAX

      
        % f_gamma_e
    
        sum_T2=0;
        sum_T4=0;
        %T2
        %sum h,j
        for k=0:1:M-2
            for j=0:1:M-k-2
                a = (p_P0*(M-1)/(p_Pj*x_k(n)))^((j+k)/2);
                b_1 = ((V+1)*p_P0*p_Pj*x_k(n)/((M-1)*(p_P0+p_Pj*x_k(n))^2))^(M/2);
                b_2 = ((V+1)*p_P0*p_Pj*x_k(n)/((M-1)*(p_P0+p_Pj*x_k(n))^2))^((M-2)/2);
                c = (1+p_Pj*x_k(n)/p_P0)^k;
                d = 2*V*sqrt((M-1)*p_P0*p_Pj*x_k(n))/(p_P0+p_Pj*x_k(n));
                e = exp(-V*(p_P0+(M-1)*p_Pj*x_k(n))/(p_P0+p_Pj*x_k(n)));
                f = pochhammer(2+k-M,j);
                sum_T2 = sum_T2 ...
                       + (-1)^j * (M-1)*p_Pj /(p_P0*factorial(j)) ...
                       * a * b_1 * c * besseli(2+j+k-M,d) * (k+1) * e * f ... //1
                       + (-1)^j * M * (M-1)* (p_P0-p_Pj*x_k(n))/(2*p_P0*x_k(n)*factorial(j)) ...
                       * a * b_1 * c * besseli(2+j+k-M,d) * e * f ... //2
                       + (-1)^(j+1) * (j+k) * (1+V) * p_Pj / (2 * p_P0 * x_k(n) * factorial(j)) ...
                       * a * b_2 * c * besseli(2+j+k-M,d) * e * f ... //3
                       + (-1)^(j+1) * V * (M-1) * (M-2) * p_Pj / ( (p_P0 + p_Pj * x_k(n)) * factorial(j) ) ...
                       * a * b_1 * c * besseli(2+j+k-M,d) * e * f ... //4
                       + (-1)^j * V * p_Pj * (M-1)^2 * (p_P0 - p_Pj * x_k(n)) ...
                       / (2 * sqrt((M-1)*p_P0 * p_Pj *x_k(n))*(p_P0 + p_Pj * x_k(n))*factorial(j)) ...
                       * a * b_1 * c * (besseli(1+j+k-M,d) + besseli(3+j+k-M,d)) * e * f;
            end
        end
        
        T2 = 2*(1-q)*q ...
           * (V*(M-1)*p_P0*p_Pj/(p_P0+p_Pj*x_k(n))^2 ...
           * exp(-V*(p_P0 + (M-1)*p_Pj*x_k(n))/(p_P0+p_Pj*x_k(n))) ...
           * hypergeom([],1,V^2*(M-1)*p_P0*p_Pj*x_k(n)/(p_P0+p_Pj*x_k(n))^2) ...
           + V*p_P0*p_Pj/(p_P0+p_Pj*x_k(n))^2 ...
           * (marcumq(sqrt(2*V*p_P0/(p_P0+p_Pj*x_k(n))),sqrt(2*V*(M-1)*p_Pj*x_k(n)/(p_P0+p_Pj*x_k(n))),2) ...
           - marcumq(sqrt(2*V*p_P0/(p_P0+p_Pj*x_k(n))),sqrt(2*V*(M-1)*p_Pj*x_k(n)/(p_P0+p_Pj*x_k(n))),1)) ...
           + sum_T2);
        %T4
        %sum k,j
        for k=0:1:M-1
            for j=0:1:M-k-1
                a = (p_P0 * (M-1)/(p_Pj * x_k(n)))^((j+k)/2);
                b_1 = (p_P0 * p_Pj * x_k(n)/((M-1)* (p_P0 + p_Pj * x_k(n))^2))^(M/2);
                b_2 = (p_P0 * p_Pj * x_k(n)/((M-1)* (p_P0 + p_Pj * x_k(n))^2))^((M-2)/2);
                c = (1 + p_Pj * x_k(n)/ p_P0)^k;
                d = 2 * V * sqrt(2 * (M-1) * p_P0 * p_Pj * x_k(n)/(p_P0 + p_Pj * x_k(n)));
                e = exp(-(M + (V+2) * p_P0) / (p_P0 + p_Pj * x_k(n)));
                f = pochhammer(1+k-M,j);
                
                sum_T4 = sum_T4 ...
                       + (-1)^j * 2^((M-2-j-k)/2) * (M-1) * sqrt(V*(M-1))*p_Pj/(p_P0 * factorial(j)) ...
                       * a * b_1 * c * besseli(2+j+k-M,d) * e * f ... //1
                       + (-1)^(j+1) * 2^((M-j-k-4)/2) * (j+k) * (M-1) * sqrt(V*(M-1)) * p_Pj ...
                       / (p_P0 * factorial(j)) ...
                       * a * b_1 * c * besseli(2+j+k-M,d) * e * f ... //2
                       + (-1)^j * 2^((M-2-j-k)/2) * sqrt(V * (M-1)) * (M + (V+2) * p_P0) * p_Pj^3 * x_k(n)^2 ...
                       / ((p_P0 + p_Pj * x_k(n))^4 * factorial(j)) ...
                       * a * b_2 * c * besseli(2+j+k-M,d) * e * f ... //3
                       + (-1)^j * 2^((M-j-k-2)/2) * k * sqrt(V * (M-1)) * p_Pj^3 * x_k(n)^2 ...
                       / ((p_P0 + p_Pj * x_k(n))^3 * factorial(j)) ... 
                       * a * b_2 * c * besseli(2+j+k-M,d) * e * f ... //4
                       + (-1)^j * 2^((M-j-k-4)/2) * sqrt(V * (M-1)) * M * p_Pj^2 * x_k(n) ...
                       / ((p_P0 + p_Pj * x_k(n))^3 * factorial(j)) ...
                       * a * b_2 * c * besseli(2+j+k-M,d) * e * f ... //5
                       + (-1)^j * 2^((M-j-k-3)/2) * e * (M-1)^2 * (V*(M-1))^(3/2) * p_Pj * f ...
                       * sqrt((p_P0 + p_Pj * x_k(n))/((M-1)*p_P0 * p_Pj * x_k(n))) ...
                       / (factorial(j)) ...
                       *(besseli(1+j+k-M,d) + besseli(3+j+k-M,d));
            end
        end
        
        T4 = q^2 * ( V^2 * (M-1)^2 * p_P0 * p_Pj^2 * x_k(n) / (p_P0 + p_Pj * x_k(n))^3 ...
           * exp(-(V*(2*p_P0 + (M-1)*p_Pj * x_k(n))/(p_P0 + p_Pj * x_k(n)))) ...
           * hypergeom([], 2, 2 * V^2 * (M-1) * p_P0 * p_Pj * x_k(n) / (p_P0 + p_Pj * x_k(n))^2) ...
           + 2 * V * p_P0 * p_Pj/(p_P0 + p_Pj * x_k(n))^2 ...
           * (marcumq(2*sqrt(V*p_P0/(p_P0 + p_Pj * x_k(n))), sqrt(2*V*(M-1)*p_Pj * x_k(n)/(p_P0 + p_Pj * x_k(n))),3) ...
           - marcumq(2*sqrt(V*p_P0/(p_P0 + p_Pj * x_k(n))), sqrt(2*V*(M-1)*p_Pj * x_k(n)/(p_P0 + p_Pj * x_k(n))),3)) ...
           + sum_T4);
    
        f_gamma_e =  T2+T4


    
        % F_gamma_d
        sum_main = 0;
        alpha = sqrt(2*V);
        beta = sqrt(2*(1+V)*p_Pj*(2^(2*R)*(1+x_k(n))-1)/(p_P0*Omega));
        c = 2*sqrt((1+V)*V*(M-1)/Omega);
        p = (1+V)/Omega;
        tilde_p = 2*p+beta^2;
        mu = M-1;
    
    %sum part
    for k=0:1:mu-1
        for j=0:1:mu-k-1
            sum_main = sum_main + (alpha/c)^(j+k) * (tilde_p/beta)^k * (-beta)^j ...
                * pochhammer(-mu+k+1 , j)/factorial(j)...
                * besseli(k+j,c*alpha*beta/tilde_p);
        end
    end
    
    
    % single channel
        F_gamma_d_k = 1 ...
                         - ((V+1)/Omega)^(M/2) * (V*(M-1))^(-(M-2)/2)*exp(-(M-1)*V) ...
                         * ((c/2)^(mu-1)*p^(-mu)*exp(c^2/(4*p))...
                         + 2/c*(c/tilde_p)^(mu)*exp((c^2/2-p*alpha^2)/tilde_p)...
                         * sum_main ...
                         - 1/p * (c/(2*p))^(mu-1)*exp(c^2/(4*p))...
                         * marcumq(beta*c/sqrt(2*p*tilde_p),alpha*sqrt(2*p/tilde_p),mu));
    % Totally
        F_gamma_d = 1 - (1- F_gamma_d_k)^M*(1- F_gamma_d_k)
        pr_secrecy = pr_secrecy + w_k(n)*F_gamma_d * f_gamma_e * exp(x_k(n))
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
        
        C_d = 1/2 * log2(1+gamma_D);
        
        C_e = 1/2 * log2(1+gamma_E);
        
        C_s = max(C_d-C_e,0);
        if (C_s < R)
            counter=counter+1;
        end
        
    end
    
    pr_secrecy_simulation(Pindex) = counter/N_max
    
end