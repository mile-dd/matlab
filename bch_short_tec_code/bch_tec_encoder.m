% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCH Encoder
function codeword = bch_tec_encoder(G, data_bits, edac)

codeword = data_bits * G;
codeword = mod(codeword, 2);

if strcmp(edac, 'dec_ted')
  % calculate an extra overall parity-check bit over all data & parity bits
  overall_parity = mod(sum(codeword), 2);
  
  % append overall parity bit to the codeword
  codeword = [overall_parity codeword];
end

end

