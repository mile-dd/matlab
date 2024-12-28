% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulation script tests BCH encoder/decoder functions
% BCH code can be DEC, DEC-TED or TEC code

clear all;
clc;

n_data_bits = 32;
num_BCH_blocks = 2000;

% Paramter 'edac' can be set to : 'tec', 'dec_ted' or 'dec';
edac = 'tec';

% BCH code structure
s_bch_edac = bch_tec_struct(n_data_bits, edac);

fprintf('Number of data bits : %d\n', n_data_bits);
fprintf('Number of data words: %d\n', num_BCH_blocks);

num_corr_err = 0;
num_uncorr_err = 0;

% Generate numbers 0 to 3 with given probabilities
num_inj_errors = randsrc(1, num_BCH_blocks, [0 1 2 3; 0.1 0.2 0.3 0.4]);

num_total_inj_err = zeros(4, 1);

for i = 1 : num_BCH_blocks
  
  % random data bits
  x = randi([0 1], n_data_bits, 1)';
  
  % BCH Encoder
  enc_bits = bch_tec_encoder(s_bch_edac.G, x, edac);
  
  % Log the number of injected errors
  num_total_inj_err(num_inj_errors(i) + 1) = num_total_inj_err(num_inj_errors(i) + 1) + 1;
  
  % Generate random bit errors
  % Single Bit Error (SBE)
  % Double Bit Error (DBE)
  % Triple Bit Error (TBE)
  random_errors = randerr(1, length(enc_bits), num_inj_errors(i));
  
  % Injection of random errors in encoded bits
  rx_bits = bitxor(enc_bits, random_errors);
  
  % BCH Decoder
  [y, corr_err, uncorr_err] = bch_tec_decoder(s_bch_edac, rx_bits, edac);
  
  % Update the number of corrected/uncorrected errors
  num_corr_err = num_corr_err + corr_err;
  num_uncorr_err = num_uncorr_err + uncorr_err;
  
  if mod(i,100) == 0
    fprintf('BCH blocks decoded : %d/%d\n', i, num_BCH_blocks);
  end
  
end

% Print reports
fprintf('\n');
fprintf('Number of corrected errors : %d \n', num_corr_err);
fprintf('Number of uncorrected errors : %d \n\n', num_uncorr_err);
fprintf('Number of error-free words : %d \n', num_total_inj_err(1));
fprintf('Number of injected SBE errors : %d \n', num_total_inj_err(2));
fprintf('Number of injected DBE errors : %d \n', num_total_inj_err(3));
fprintf('Number of injected TBE errors : %d \n', num_total_inj_err(4));

fprintf('\n');

if strcmp(edac, 'tec')
  if sum(num_total_inj_err) == num_BCH_blocks
    fprintf('Success. All injected errors have been corrected !\n');
  end
end

if strcmp(edac, 'dec_ted')
  if sum(num_total_inj_err(1:3)) == num_corr_err
    fprintf('Success. All injected errors (SBE + DBE) have been corrected !\n');
  end
  if (num_total_inj_err(4) == num_uncorr_err)
    fprintf('Success. All TBEs have been correctly detected !\n');
  end
end

if strcmp(edac, 'dec')
  if sum(num_total_inj_err(1:3)) == num_corr_err
    fprintf('Success. All injected errors (SBE + DBE) have been corrected !\n');
  end
  if (num_total_inj_err(4) ~= num_uncorr_err)
    warning('Some TBEs have been not correctly detected - expected behavior for DEC code');
    fprintf('Those errors are misinterpreted as single/double bit errors\n');
  end
end


