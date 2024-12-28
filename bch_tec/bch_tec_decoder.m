% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s_dec = bch_tec_decoder(received_code, s_bch)

% Three-Error-Correcting (TEC) BCH Decoder, t = 3;
% Double-Error-Correction and Single-Error-Correcting codes 
% are also supported, (t = 1, t = 2);
% 
% It uses Peterson algoritm to calculate error location polynomial
% instead of using BM (Berlekamp邦assey) algorithm;
% For the root computation of error locator polynomial, instead of
% Chien search algorithm, it uses another method
% based on LUT (lookup table)
%
% it processes one BCH codeword
% received_code: received bits
% s_bch: BCH code structure, see bch_code_struct.m

received_code = uint32(received_code);

s = s_bch;
clear s_bch;

n = s.n;  % codeword length in bits
t = s.t;  % error correction capability
k = s.k;  % data length in bits
%=================================================
% 1. syndrome calculation
syndrome = uint32(zeros(1, 2*t));

for i = 1:n
  syndrome = bitxor(gf_mul(syndrome, s.alpha_synd, s), repmat(received_code(i),1,2*t));
end
%=================================================
% 2. error locator polynomial (lambda) calculation
lambda = peterson_algorithm(syndrome, s);
index = find(lambda(2:t+1));

% find lambda degree
if isempty(index)
  lambda_degree = 0;
else
  lambda_degree = index(length(index));
end
%=================================================
% 3. Find error locations
error_location = roots_calculation(lambda, s);
number_roots = sum(error_location > 0);
%=================================================
% 4. error correction in received code
full_corrected_code = received_code;
full_corrected_code(error_location(error_location > 0)) = ...
  1 - received_code(error_location(error_location > 0));

if lambda_degree == number_roots
  % correct received symbols
  corrected_code = full_corrected_code(1:k);
else % deliver original received code
  corrected_code = received_code(1:k);
end

% save in structure s_dec
s_dec.corrected_code = corrected_code;
s_dec.syndrome = syndrome;
s_dec.lambda = lambda;
s_dec.error_location = error_location;
s_dec.lambda_degree = lambda_degree;
s_dec.number_roots = number_roots;

