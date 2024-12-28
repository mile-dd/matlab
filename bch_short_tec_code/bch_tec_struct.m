% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCH code structure, it contains :
% Correctable error patterns table 
% Syndrome table for all correctable error patterns
% Generator matrix in systematic form
% Parity-check matrix in systematic form

function s = bch_tec_struct(n_data_bits, edac)

% full parity-check matrix
[H] = generate_bch_matrix_cyclic(n_data_bits, edac);

if strcmp(edac, 'dec')
  t = 2;
elseif strcmp(edac, 'dec_ted')
  t = 2;
elseif strcmp(edac, 'tec')
  t = 3;
end

% size of H matrix
[num_row, num_col] = size(H);

% H_ted: This matrix provides DEC-TED (Triple Error Detection)
% Hamming distance is increased for 1, so d = 6;
% To get H_ted, extend H in the following way:
%
% H_ted = |1 1 1 . . . 1|
%         |0            |
%         |.            |
%         |.     H      |
%         |.            |
%         |0            |
H_ted = zeros(num_row + 1, num_col + 1);
for i = 1: num_row
  for j = 1:num_col
    H_ted(i+1,j+1) = H(i,j);
  end
end
% first row shall be all '1' vector
H_ted(1, :) = ones(1, num_col + 1);

% number of data(source) bits
K = num_col - num_row;

% total number of bits: data + parity bits
N = num_col;

H_parity = H(:, num_row+1:N);

% Generator matrix in the systematic form
G = cat(2, H_parity', eye(K));

if strcmp(edac, 'dec_ted')
  % For DEC-TED code, total number of bits is greater for one;
  N = N + 1;
end

% Number of correctable error patterns
if t == 2
  n_err_pattern = N*(N+1)/2;
elseif t == 3
  n_err_pattern = N*(N+1)/2 + N*(N-1)*(N-2)/6;
end

% Initialize error pattern table
err_pattern_table = zeros(n_err_pattern, N);

% error patterns for single bit errors: N
err_pattern_table(1:N,1:N) = eye(N, N);

% error patterns for double bit errors: N*(N-1)/2
new_entry = N + 1;
for m = 1 : N-1
  for n =  m+1 : N
     err_pattern_table(new_entry, m) = 1;
     err_pattern_table(new_entry, n) = 1;
     new_entry = new_entry + 1;
  end
end

% error patterns for triple bit errors: N*(N-1)*(N-2)/6;
if t == 3  
  new_entry = N*(N+1)/2 + 1;
  for m = 1 : N-2
    for n = m+1 : N-1
      for i = n+1 : N	
        err_pattern_table(new_entry, m) = 1;
        err_pattern_table(new_entry, n) = 1;
        err_pattern_table(new_entry, i) = 1;
        new_entry = new_entry + 1;
      end
    end
  end
end

% Calculate the syndromes for all correctable error patterns
if strcmp(edac, 'dec_ted')
  syndrome_table = double(err_pattern_table) * double(H_ted');
else
  syndrome_table = double(err_pattern_table) * double(H');
end

syndrome_table = mod(syndrome_table, 2);
% % convert to decimal integer
syndrome_table_de =  bi2de(syndrome_table, 'left-msb');

% save in structure
s.H = double(H);
s.H_ted = double(H_ted);
s.G = double(G);
s.err_pattern_table = double(err_pattern_table);
s.syndrome_table_de = syndrome_table_de;

end

