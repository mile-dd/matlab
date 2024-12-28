% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [error_location] = chien_search_full(lambda, s_bch)

% Chien search algorithm to find error positions in the received code
% Code length is always set to n_max to have complete search !

% lambda: error location polynomial
% s_bch: BCH code structure, see bch_code_struct.m

s = s_bch;
clear s_bch;

t = s.t; % error correction capability
n = s.n_max; % code length is always set to n_max to have complete search !
n_max = s.n_max; % maximum code length. It must be of form 2^m-1

% accelerator lookup tables using function handle
apo = @(x) s.alpha_power(1+uint32(x));

alpha_cs = apo(0:t);

% Initialization of Chien cells
% Correct initialization is needed for a shortened version of BCH code
n_index = n_max - n;
cs_cells(1:t+1) = gf_mul(lambda(1:t+1), apo(mod(n_index*(0:t), n_max)), s);

index = 1;
zero = uint32(0);
error_location = zeros(1, t);

for ii = 1:n
  cs_cells = gf_mul(cs_cells, alpha_cs, s);
  lambda_value = zero;
  
  for jj = 1:t+1
    lambda_value = bitxor(lambda_value, cs_cells(jj));
  end
  
  if (lambda_value == 0) % an error location has been found
    error_location(index) = ii;
    index = index + 1;
  end
  
end

end

