% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulation script tests BCH encoder/decoder functions
% To set BCH code parameters, it can be used Table 6.1 (page 195)
% in book "Error Control Coding" (Authors : Shu Li, Daniel J. Costello), Second Edition, 2004.
% It can be also used MATLAB function "bchnumerr(N)"

clear all;
clc;

% Parameters:
n_max = 1023; % maximum number of bits per BCH codeword. 
% n_max must be of form 2^m-1;
m = log2(n_max+1);
n = 768; % actual number of bits per BCH codeword (code length)
prim_poly = primpoly(m); % primitive polynomial for GF(2^m) = primpoly(m)

bch_param = bchnumerr(n_max);
[n_bch_param, ~] = size(bch_param); 

% n_param_set belongs to interval [1 n_bch_param];
% for shortened BCH code (n < n_max), n_param_set must be set to some value 
% less than n_bch_param to provide that number of data bits 
% per BCH codeword is greater than zero !
n_param_set = 40; 

num_BCH_blocks = 1000;
num_BCH_conf = 20;

for ii = 1 : num_BCH_conf

  n_random = randsrc(1,1,[1:n_param_set]);
  % take k_max from "bch_param" matrix   
  % Values n_max and k_max must produce a valid narrow-sense BCH
  k_max = bch_param(n_random, 2);
  
  % Create BCH code structure 
  s_bch = bch_code_struct(n_max, k_max, n, prim_poly);
  k = s_bch.k; % number of data bits per BCH codeword
  t = s_bch.t; % error correction capabilty
  
  fprintf('global iteration %d/%d  \n', ii, num_BCH_conf);  
  fprintf('error correction capability : %d bits \n', t);

  check = zeros(1,num_BCH_blocks);

  for i = 1 : num_BCH_blocks

    % random data bits
    x = randi([0 1], k, 1)';

    % BCH Encoder
    enc_bits = bch_encoder(x, s_bch);
    
    n_random_errors = randsrc(1,1,[1:t]);    
    random_errors = randerr(1, n, n_random_errors); 
    % Injection of random errors in encoded bits
    received_code = bitxor(enc_bits, uint32(random_errors));

    % BCH Decoder
    s_dec = bch_decoder(received_code, s_bch);

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

end
