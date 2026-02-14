%% Simulacao simples de OFDM no dominio do tempo
%% Autor: Francisco Muller

close all;
clear all;
%% Parametros de entrada
pkg load communications % necessario caso Octave

N = 64; %Tamanho FFT - quantidade de subcanais. Nao pode ser menor que n_cp
M = 16; %Ordem da modulacao (QAM)
N_ofdm = 10000; % number of OFDM symbols
%n_cp = 0; % tamanho do prefixo ciclico. Idealmente igual ou maior que length(h)-1
f_c = 10^5; % frequencia da portadora em Hz (apenas para plotar resposta do canal)

lcp = input("Usará a largura do canal como número de prexixo cíclico? [S] | N : ", "s");

% Escolha a SNR media aqui
%SNR_time_domain = 20; %(dB)

SNR_time_domain = input("Digite o valor inteiro da SNR (dB) ", "s");
SNR_time_domain = str2double(SNR_time_domain);

if (isnan(SNR_time_domain) || SNR_time_domain <= 0)
    disp("Entrada irregular, aplicando SNR padrão de 20dB.");
    SNR_time_domain = 20;
endif

disp("SNR: ");
disp(SNR_time_domain);


h = [-1.4190 - 0.2916i
        0 +      0i
        0 +      0i
   0.2980 + 0.5857i
        0 +      0i
  -0.0212 - 0.0836i
  -0.0246 - 0.0758i
        0 +      0i
  -0.0094 - 0.0554i];

if (lcp == 'n'|| lcp == 'N' )
  n_cp = 0;
else
  n_cp = length(h);
end


% --- BLOCO DE CRIAÇÃO DE PASTA E NOME ---
% Define o nome da pasta
dir_name = sprintf('simulation_SNR%d_NCP%d', SNR_time_domain, n_cp);

% Cria a pasta se ela não existir
if ~exist(dir_name, 'dir')
    mkdir(dir_name);
end

% Define um prefixo padrão para os arquivos para economizar digitação
file_prefix = fullfile(dir_name, sprintf('SNR_%d_NCP%d_', SNR_time_domain, n_cp));

H = fft([h; zeros(N-length(h),1)]); % funcao de transferencia
%Plot do canal
dt = 10^-3;
stem([0:(length(h)-1)]*dt,h);
title('Resposta ao Impulso');
xlabel('Time (s)');
ylabel('h(t)');

%para salvar o gráfico automaticamente
print(sprintf('%sResposta_Impulso.png', file_prefix), '-dpng');

df = 1 / (N * dt);
f = [ -(ceil((N-1)/2):-1:1), 0, (1:floor((N-1)/2)) ] * df + f_c;
figure, plot(f,10*log10(abs(fftshift((H)))));
title('Magnitude da Funcao de Transferencia');
xlabel('Frequency (Hz)');
ylabel('{|H(f)}| (dB)');

print(sprintf('%sMagnitude_Transferencia.png', file_prefix), '-dpng');


figure, plot(f,arg(fftshift((H))));
title('Fase da Funcao de Transferencia');
xlabel('Frequency (Hz)');
ylabel('\angle{H(f)}');

print(sprintf('%sFase_Transferencia.png', file_prefix), '-dpng');


%% Simbolos de entrada
x = randint(N,N_ofdm,M); %Gerando numeros de 0 a M-1
s = qammod(x,M); %Gerando simbolos M-QAM
E_s = meansq(s(:)); % Energia media dos simbolos
s_norm = s/sqrt(E_s); % simbolos normalizados (E_media=1)
s_ifft = sqrt(N)*ifft(s_norm); %passando pela IDFT normalizada

% Adicionando prefixo ciclico
s_cp = s_ifft(N-n_cp+1:N,:);
%s_zp = zeros(n_cp,N_ofdm);
s_with_cp = [s_cp; s_ifft]; % adicionando prefixo ciclico.
%s_with_zp = [s_zp; s_ifft]; % adicionando zero padding.

% Preparando o sinal para a convolucao com o canal
s_serial = s_with_cp(:); % matriz vira linha
%s_serial_zp = s_with_zp(:);

% plot sinal transmitido (1000 amostras)
    %Sinal Transmitido Real
figure, plot([1:1000]*dt,real(s_serial(1:1000)));
    % Título adicionado
title('Sinal Transmitido - Real');
print(sprintf('%sSinal_Tx_Real.png', file_prefix), '-dpng');
    %Sinal Transmitido Imaginário
figure, plot([1:1000]*dt,imag(s_serial(1:1000)));
title('Sinal Transmitido - Imag');
print(sprintf('%sSinal_Tx_Imag.png', file_prefix), '-dpng');

figure, pwelch(s_serial, 'shift','dB');
title('PSD Sinal Serial');
print(sprintf('%sPSD_Tx.png', file_prefix), '-dpng');
%figure, pwelch(s_serial_zp, 'shift','dB');


%% Passando pelo canal
s_conv = conv(s_serial,h);
figure, pwelch(s_conv, 'shift','dB');
title('PSD Sinal Recebido');
print(sprintf('%sPSD_Rx.png', file_prefix), '-dpng');

%s_conv_zp = conv(s_serial_zp,h);
%figure, pwelch(s_conv_zp, 'shift','dB');
% Adicionar AWGN (com base no sinal no tempo)
s_ruido = awgn(s_conv,SNR_time_domain,'measured','dB');
% TOD: Adicionar PSD ruido

%% Recepcao
r_tmp = s_ruido(1:end-length(h)+1); % tirando o "resto" da convolucao
r_with_cp = reshape(r_tmp,n_cp+N,N_ofdm);
r = r_with_cp(n_cp+1:end,:); % retirando prefixo
z = fft(r)/sqrt(N); % Aplicando DFT normalizada
z_norm = z*sqrt(E_s);

% one-tap equalization
H_eq = 1./repmat(H,1,N_ofdm);
z_eq = (H_eq.*z_norm); %equalizacao

%deteccao
x_est = qamdemod(z_eq,M);

% Total de bits transmitidos
NUM_BIT_Total = log2(M)*N*N_ofdm

%% Bit and Symbol Error Rate Calculation
[NUM_BIT, BER] = biterr(x, x_est)
[NUM_SYM, SER] = symerr(x, x_est)

BER_sub = zeros(1,N);
NUM_BIT_sub = zeros(1,N);
SER_sub = zeros(1,N);
NUM_SYM_sub = zeros(1,N);

for nn = 1:N
  [NUM_BIT_sub(nn), BER_sub(nn)] = biterr(x(nn,:),x_est(nn,:));
  [NUM_SYM_sub(nn), SER_sub(nn)] = symerr(x(nn,:),x_est(nn,:));
end



% Plotar (ou nao) constelacoes

% Melhor e pior subcanais (em termos de magnitude)
[sub_best, sub_best_idx]=max(abs(H))
[sub_worst, sub_worst_idx]=min(abs(H))

%scatterplot para o melhor subcanal
z_plot = scatterplot(z_eq(sub_best_idx,:),1,0,"r+");
hold on, grid on;
scatterplot(s(sub_best_idx,:),1,0,"b+",z_plot);
title('Constelacao QAM - Melhor Subcanal');
print(sprintf('%sMelhor_Canal.png', file_prefix), '-dpng');

%scatterplot para o melhor subcanal sem equalizacao
z_plot = scatterplot(z_norm(sub_best_idx,:),1,0,"r+");
hold on, grid on;
scatterplot(s(sub_best_idx,:),1,0,"b+",z_plot);
title('Constelacao QAM - Melhor Subcanal sem EQ');
print(sprintf('%sMelhor_Canal_SemEQ.png', file_prefix), '-dpng');

%scatterplot para o pior subcanal
z_plot = scatterplot(z_eq(sub_worst_idx,:),1,0,"r+");
hold on, grid on;
scatterplot(s(sub_worst_idx,:),1,0,"b+",z_plot);
title('Constelacao QAM - Pior Subcanal');
print(sprintf('%sPior_Canal.png', file_prefix), '-dpng');

%scatterplot para o pior subcanal sem equalizacao
z_plot = scatterplot(z_norm(sub_worst_idx,:),1,0,"r+");
hold on, grid on;
scatterplot(s(sub_worst_idx,:),1,0,"b+",z_plot);
title('Constelacao QAM - Pior Subcanal sem EQ');
print(sprintf('%sPior_Canal_SemEQ.png', file_prefix), '-dpng');

%scatterplot para o subcanal 10 sem equalizacao
z_plot = scatterplot(z_norm(10,:),1,0,"r+");
hold on, grid on;
scatterplot(s(10,:),1,0,"b+",z_plot);
title('Constelacao QAM - Subcanal 10 sem EQ');
print(sprintf('%sCanal_10_SemEQ', file_prefix), '-dpng');

%scatterplot para subcanal 10
z_plot = scatterplot(z_eq(10,:),1,0,"r+");
hold on, grid on;
scatterplot(s(10,:),1,0,"b+",z_plot);
title('Constelacao QAM - Subcanal 10');
print(sprintf('%sCanal_10.png', file_prefix), '-dpng');


%% Bit and Symbol Error Rate Calculation (pior subcanal)
[NUM_BIT_worst, BER_worst] = biterr(x(sub_worst_idx,:), x_est(sub_worst_idx,:))
[NUM_SYM_worst, SER_worst] = symerr(x(sub_worst_idx,:), x_est(sub_worst_idx,:))

%% Bit and Symbol Error Rate Calculation (melhor subcanal)
[NUM_BIT_best, BER_best] = biterr(x(sub_best_idx,:), x_est(sub_best_idx,:))
[NUM_SYM_best, SER_best] = symerr(x(sub_best_idx,:), x_est(sub_best_idx,:))

%% Bit and Symbol Error Rate Calculation (subcanal 10)
[NUM_BIT_10, BER_10] = biterr(x(10,:), x_est(10,:))
[NUM_SYM_10, SER_10] = symerr(x(10,:), x_est(10,:))

