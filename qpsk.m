% Constellations and we plot the resulting capacity against SNR
clear
clc
%QPSK
h_yx = 0.5*log2(2*pi*exp(1));
% For calculating h(Y) we generate the bits i.i.d with p = 1/2 and compute
% the capacity numerically using the law of large numbers
snr_db = (-3:20)';
snr = 10.^(snr_db/10);
num_bits = 1e6;
cap_qpsk = zeros(length(snr),1);
for ctr = 1:length(snr)
    bits = round(rand(num_bits,1)); % Generating the bits with equal probability
    bpsk_symbols = 1 - 2*bits; % Mapping bits to symbols
    A = sqrt(snr(ctr));
    noise = randn(num_bits,1);
    y = A*bpsk_symbols + noise;
    prob_y = (0.5/sqrt(2*pi))*(exp(-0.5*(y - A).^2) + exp(-0.5*(y + A).^2));
    h_y = -sum(log2(prob_y))/num_bits;
    % Capacity of QPSK = 2*Capacity of BPSK
    cap_qpsk(ctr) = 2*(h_y - h_yx);
end

% Plotting the capacity as a function of SNR
plot(snr_db,cap_qpsk,'--r')
hold on
ylim([0 5])

% ---------------------------------------------------------------------- %
% ---------------- CAPACITY FOR 16 - QAM ------------------------------- %

cap_qam = zeros(length(snr),1);
num_symbols = 1e6;
for ctr = 1:length(snr)
    d = sqrt(snr(ctr)/5);
    % Generating the symbols
    seq_1_3 = 2*round(rand(num_symbols,1)) + 1;
    seq_1_neg1 = 1 - 2*round(rand(num_symbols,1));
    pam_symbols = seq_1_3.*seq_1_neg1.*d;
    noise = randn(num_symbols,1);
    y = pam_symbols + noise;
    % Now we calculate h(Y) as - E[log2(p(y)]
    prob_y = (0.25/sqrt(2*pi))*(exp(-0.5*(y - d).^2) + exp(-0.5*(y + d).^2) + exp(-0.5*(y + 3*d).^2) + exp(-0.5*(y - 3*d).^2));
    h_y = -sum(log2(prob_y))/num_symbols;
    %Capacity of 16-QAM = 2*Capacity of 4-PAM
    cap_qam(ctr) = 2*(h_y-h_yx);
end
plot(snr_db,cap_qam,'-xb')
% ---------------------------------------------------------------------- %
% ------------------- CAPACITY FOR 16-PSK ------------------------------ %

h_yx = log2(pi*exp(1));
num_symbols = 1e5;
prob_y = zeros(num_symbols,1);
cap_psk = zeros(length(snr),1);
for ctr = 1:length(snr)
    A = sqrt(snr(ctr));
    levels = round(16*rand(num_symbols,1));
    psk_symbols = A*exp(j*(2*pi/16)*levels);
    noise = (1/sqrt(2))*(randn(num_symbols,1) + j*randn(num_symbols,1));
    rec_symbols = psk_symbols + noise;
    const_points = A*exp(j*(2*pi/16)*(0:15));
    const_points = const_points(:);
    rec_symbols_repeated = repmat(rec_symbols.',16,1);
    psk_symbols_repeated = repmat(const_points,1,num_symbols);
    diff_matrix = (rec_symbols_repeated - psk_symbols_repeated);
    prob_y = sum(exp(-(abs(diff_matrix)).^2),1)/(pi*16);
    h_y = -sum(log2(prob_y))/num_symbols;
    cap_psk(ctr) = (h_y - h_yx);
end
plot(snr_db,cap_psk,'-+k')
grid on
xlabel('SNR(in dB)')
ylabel('Capacity (in bits/channel use)')
legend('QPSK','16-QAM','16-PSK')
% ---------------- Plots in terms of Eb/N0 -------------------

ebno_qpsk_db = 10*log10(snr./cap_qpsk);
ebno_qam_db = 10*log10(snr./cap_qam);
ebno_psk_db = 10*log10(snr./cap_psk);
figure
plot(ebno_qpsk_db,cap_qpsk,'-*r')
hold on

plot(ebno_qam_db,cap_qam,'-xb')
plot(ebno_psk_db,cap_psk,'-+k')
grid on
xlabel('EB/N0(in dB)');
ylabel('Capacity (bits/channel use)');
axis tight
legend('QPSK','16-QAM','16-PSK')
% -------------------------------------------------------------------------