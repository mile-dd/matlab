% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulation script estimates the BER performance of BCH decoder;
% BPSK modulation over AWGN channel 
clear all; 
clc;
% Parameters:
n_max = 1023; % maximum number of bits per BCH codeword.
k_max = 863; 
% n_max must be of form 2^m-1;
m = log2(n_max+1);
n = 768; % actual number of bits per BCH codeword (code length)
prim_poly = primpoly(m); % primitive polynomial for GF(2^m) = primpoly(m)

% Create BCH code structure 
% Values n_max and k_max must produce a valid narrow-sense BCH
s_bch = bch_code_struct(n_max, k_max, n, prim_poly);
k = s_bch.k; % number of data bits per BCH codeword
t = s_bch.t; % error correction capabilty
               
EbN0 = [5.0 5.5 6.0 6.5];
L = length(EbN0);

EbN0_coded = EbN0 + 10*log10(k/n);
N0_coded = 1./(exp(EbN0_coded*log(10)/10));
N0_uncoded = 1./(exp(EbN0*log(10)/10));

% Number of BCH blocks
num_bch_blocks = 1000;

ber_bch = zeros(1, L); % BER with BCH coding 
ber_uncoded = zeros(1, L); % BER uncoded

% Output structure s_info
s_info = [ ];

fprintf('error correction capability : %d bits \n', t);

for kk = 1:length(EbN0)
  
  fprintf('global iteration %d/%d for EbN0_dB = %.2f\n\n', kk, L, EbN0(kk));
  err_bits = 0;
  err_bits_unc = 0;
  
  for nn = 1:num_bch_blocks
    
    % random data bits
    x = randi([0 1], k, 1)';
    
    % BCH Encoder
    encoded_bits = bch_encoder(x, s_bch);
       
    % BPSK modulation
    bpsk_mod = 2*double(encoded_bits) - 1;
    
    % AWGN Channel
    noise_coded = sqrt(N0_coded(kk)/2)*randn(size(bpsk_mod));
    noise_uncoded = sqrt(N0_uncoded(kk)/2)*randn(size(bpsk_mod));
    
    % Add noise
    y_coded = bpsk_mod + noise_coded;
    y_uncoded = bpsk_mod + noise_uncoded;
    
    % Hard decision
    y_coded_hard = (y_coded > 0);
    y_uncoded_hard = (y_uncoded > 0);
    
    % BCH Decoder
    s_dec = bch_decoder(y_coded_hard, s_bch);    
   
    err_bits_new = sum(s_dec.corrected_code ~= x);
    err_bits = err_bits + err_bits_new;
    
    err_bits_unc_new = sum(y_uncoded_hard(1:k) ~= x);
    err_bits_unc = err_bits_unc + err_bits_unc_new;
    
    if mod(nn, 100) == 0
      fprintf('Number of BCH blocks decoded : %d/%d\n', nn, num_bch_blocks);
    end
    
  end
  
  % BER with BCH coding
  ber_bch(kk) = err_bits/(num_bch_blocks*k);
  
  % Uncoded BER
  ber_uncoded(kk) = err_bits_unc/(num_bch_blocks*k);
  
end

% save in structure simulations results and parameters
s_info.ber_bch = ber_bch;
s_info.ber_uncoded = ber_uncoded;
s_info.EbN0 = EbN0;
s_info.s_bch = s_bch;
fprintf('Simulation has been finished, see structure s_info\n');