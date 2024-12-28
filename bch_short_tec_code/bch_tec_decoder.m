% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, corr_err, uncorr_err] = bch_tec_decoder(s_bch_tec, rx_bits, edac)

% BCH Decoder: Triple-Bit-Error correction capability (TEC)
% s_bch_dec: BCH code structure, see bch_tec_struct.m
% y: decoded data bits
% corr_err : flag indicating that error is detected and corrected
% uncorr_err : flag indicating that error is detected,
% but it is uncorrectable

if strcmp(edac, 'dec_ted')
  H = s_bch_tec.H_ted;
else
  H = s_bch_tec.H;
end

[num_row, ~] = size(H);

err_pattern_table = s_bch_tec.err_pattern_table;
syndrome_table_de = s_bch_tec.syndrome_table_de;

% calculate syndromes
syndrome = rx_bits * H';
syndrome = mod(syndrome, 2);

if (sum(syndrome) == 0)
  corr_err = 1;
  uncorr_err = 0;
  y = rx_bits(num_row+1:end);
end

if (sum(syndrome) ~= 0)
  % convert to decimal integer
  syndrome_de = bi2de(syndrome, 'left-msb');
  
  % Map the syndromes to the error pattern
  index_synd = find(syndrome_de == syndrome_table_de);
  
  if ~isempty(index_synd)
    % correct errors
    yc = bitxor(rx_bits, err_pattern_table(index_synd,:));
    % discard parity bits
    y = yc(num_row+1:end);
    corr_err = 1;
    uncorr_err = 0;
  else % index_synd is empty matrix
    corr_err = 0;
    uncorr_err = 1;
    y = rx_bits(num_row+1:end);
  end
end

end

