% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = bch_tec_code_struct(n_max, k_max, n, prim_poly)

% n_max : maximum number of bits per BCH codeword. It must be of form 2^m-1
% k_max : maximum number of data bits per BCH codeword
% Values n_max and k_max must produce a valid narrow-sense BCH
% n : actual code length in bits
% prim_poly : primitive polynomial for GF(2^m) = primpoly(m)

s = [];

m = log2(n_max+1); % integer between 3 and 16 to define GF(2^m)

n_max = 2^m - 1; % maximum number of bits per BCH codeword
k = n - (n_max - k_max); % actual number of data bits per BCH codeword

% Generator polynomial
[gen_poly, t] = bchgenpoly(n_max, k_max, prim_poly);

s.m = m;
s.n = n;
s.k = k;
s.t = t; % error-correction capabilty
s.n_max = n_max;
s.k_max = k_max;

s.prim_poly = prim_poly;
gp = uint32(gen_poly.x);
gp = fliplr(gp);
s.gp = gp;

% calculate alpha power
alpha = gf(2, m, prim_poly); % primitive element 'alpha' of GF(2^m)
ap = alpha.^(0:n_max);
alpha_power = uint32(ap.x);
s.alpha_power = alpha_power;

index_gf_element = zeros(1, n_max+1);

% find GF element in alpha_power array
for nn = 2:n_max + 1
  index_gf_element(nn) = find(nn-1 == alpha_power(1:n_max));
end

s.index_gf_element = index_gf_element;

% create alpha array
nn = 1:2*t;
alpha_synd = alpha.^nn;
alpha_synd_num = uint32(alpha_synd.x);
s.alpha_synd = alpha_synd_num;

% calculate 1/x
inv_x = uint32(zeros(1,n_max));
for kk = 1:n_max
  index_kk = find(kk == alpha_power(1:n_max));
  inv_x(kk) = alpha_power(n_max+2-index_kk);
end

square_root_x = uint32(zeros(1,n_max));
% calculate square root
for kk = 1:n_max
  index_kk = find(kk == alpha_power(1:n_max));
  square_root_x(kk) = alpha_power(1 + mod((index_kk-1)*2^(m-1), n_max));
end

% Save calculated 1/x and square root;
s.inv_x = inv_x;
s.square_root_x = square_root_x;

roots_three = zeros(n_max, t);
% pre-calculate roots of lambda = z^3 + z + c = 0,
% for each possible value of c;
% Chien search is used;
for c = 1:n_max
  lambda = [c 1 0 1];
  error_location = chien_search_full(lambda, s);
  for i = 1:t
    if error_location(i) ~= 0
      roots_three(c, i) = alpha_power(error_location(i)+1);
    end
  end
end

roots_three_exception = zeros(n_max, t);
% pre-calculate roots of lambda = z^3 + b = 0,
% for each possible value of c;
% Chien search is used;
for b = 1:n_max
  lambda = [b 0 0 1];
  error_location = chien_search_full(lambda, s);
  for i = 1:t
    if error_location(i) ~= 0
      roots_three_exception(b, i) = alpha_power(error_location(i)+1);
    end
  end
end

roots_two = zeros(n_max, 2);

for b = 1:n_max
  lambda = [b 1 1];
  error_location = chien_search_full_two(lambda, s);
  for i = 1:2
    if error_location(i) ~= 0
      roots_two(b, i) = alpha_power(error_location(i)+1);
    end
  end
end

s.roots_three = roots_three;
s.roots_three_exception = roots_three_exception;
s.roots_two = roots_two;
