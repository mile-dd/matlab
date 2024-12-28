% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulation script tests BCH encoder/decoder functions for TEC

clear all;
clc;

% Parameters:
n_max = 1023; % maximum number of bits per BCH codeword. It must be of form 2^m-1
k_max = 993;  % maximum number of data bits per BCH codeword
% Values n_max and k_max must produce a valid narrow-sense BCH
% Parameters can be set by using MATLAB function "bchnumerr(N)"

m = log2(n_max+1);
n = 1023; % actual number of bits per BCH codeword (code length)
prim_poly = primpoly(m); % primitive polynomial for GF(2^m) = primpoly(m)

% BCH code structure
s_bch = bch_tec_code_struct(n_max, k_max, n, prim_poly);
k = s_bch.k; % number of data bits per BCH codeword
t = s_bch.t; % error correction capabilty

fprintf('error correction capability : %d bits \n', t);

num_BCH_blocks = 1000;

check = zeros(1,num_BCH_blocks);

% Generate numbers 0 to 3 with given probabilities
num_err = randsrc(1, num_BCH_blocks, [0 1 2 3; 0.1 0.2 0.3 0.4]);

for i = 1 : num_BCH_blocks
  % random data bits
  x = randi([0 1], k, 1)';
  
  % BCH Encoder
  enc_bits = bch_encoder(x, s_bch);
  
  random_errors = randerr(1, n, num_err(i));
  % Injection of random errors in encoded bits
  received_code = bitxor(enc_bits, uint32(random_errors));

  % BCH Decoder, Three-Error-Correcting
  s_dec = bch_tec_decoder(received_code, s_bch);
  
  if mod(i,100) == 0
    fprintf('BCH blocks decoded : %d/%d\n', i, num_BCH_blocks);
  end
  
  check(i) = isequal(x, s_dec.corrected_code);
  
  if check(i) ~= 1
    error('BCH Decoder: Errors are not corrected');
  end
end

if sum(check) == num_BCH_blocks
  fprintf('Simulation has been finished. No Errors !\n');
end

