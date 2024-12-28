% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enc_bits = bch_encoder(x, s_bch)

% BCH Encoder
% x : data bits
% enc_bits : encoded bits = [x parity_bits]
% s_bch : BCH code structure, see bch_code_struct.m

s = s_bch;
clear s_bch;

% Generator polynomial
gp = s.gp;
L = length(gp);

% Number of data bits per BCH codeword
k = s.k;

z = uint32(zeros(1,L-1));
for i_bit = 1:k
  % GF Adder
  xx = bitxor(z(L-1), x(i_bit));
  % multiplication with the generator polynominal coefficients  
  v = xx.*gp(1:L-1);
  % GF Adders
  z = bitxor(v, [uint32(0) z(1:L-2)]);
end

parity_bits = fliplr(z);

enc_bits = [x parity_bits];

end

