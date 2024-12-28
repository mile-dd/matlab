% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lambda = RiBM_BCH_algorithm(syndrome, s_bch)

% calculation of error location polynomial (lambda) using RiBM algorithm. 
% RiBM algorithm is optimized for BCH code. 
% The number of iterations is 't' which is two times less comparing
% to the RiBM algorithm for RS code.

% Reference literature:
% Xinmiao Zhang, "VLSI Architecture for Modern Error-Correcting Codes", 2016.
% page 209-211, Algorithm 21, Figure 7.4

% s_bch: BCH code structure, see bch_code_struct.m

s = s_bch;
clear s_bch;

t = s.t;

% Very first initialization of "delta" and "theta" variables
delta = uint32([zeros(1, 3*t),   1, 0, 0]);
theta = uint32([zeros(1, 3*t-1), 1, 0]);

one = uint32(1);
zero = uint32(0);

gamma = one;
kk = -1;

% Initialize with syndromes
theta(1:2*t) = [syndrome(2:2*t) zero];
delta(1:2*t) = syndrome;

for r = 1:t 
   
  % step RiBM.1:
  delta_0 = delta;
  delta(1:3*t+1) = bitxor(gf_mul(gamma, delta_0(3:3*t+3), s), gf_mul(delta_0(1), theta, s));
  
  % step RiBM.2:
  if((delta_0(1) ~= zero) && (kk >= -1))
    theta(1:3*t) = delta_0(3:3*t+2);
    theta(3*t+1) = zero;
    gamma = delta_0(1);
    kk = -kk-2;
  else
    kk = kk+2;
  end
end
% delta
lambda = delta(t+1:2*t+1);  % error locator polynominal

end

