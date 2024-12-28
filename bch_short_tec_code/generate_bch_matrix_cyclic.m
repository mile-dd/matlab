% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Parity-Check Matrix of BCH code in systematic form
% Cyclic property of BCH code is used for the matrix generation

function H = generate_bch_matrix_cyclic(n_data_bits, edac)

  switch edac
    case 'dec'
      fprintf('BCH code is Double-Error-Correction (DEC)\n');
      t = 2;
    case 'dec_ted'
      fprintf('BCH code is Double-Error-Correction, Triple-Error-Detection (DEC-TED)\n');
      t = 2;  
    case 'tec'
      fprintf('BCH code is Triple-Error-Correction (TEC)\n');
      t = 3;
    otherwise
      error('Error: EDAC option is unknown');
  end

  switch n_data_bits
    % set BCH parameters (n, k)
    case 16
      n = 31;
      if t == 2
        k = 21;
      elseif t == 3
        k = 16; 
      end
      m = log2(n+1);
      n_delete = k - n_data_bits;
      
    case 32
      n = 63;
      if t == 2
        k = 51;
      elseif t == 3
        k = 45; 
      end
      m = log2(n+1);
      n_delete = k - n_data_bits;
      
    case 64 
      n = 127;
      if t == 2
        k = 113;
      elseif t == 3
        k = 106;
      end
      m = log2(n+1);
      n_delete = k - n_data_bits;
      
    otherwise
      error('Error : number of data bits not supported');
  end
  % generator polynomial of BCH for parameters n, k, m
  % error correction capability is t = 2 or t = 3;
  [gen_poly, ~] = bchgenpoly(n, k, primpoly(m));
  
  fprintf('\n');
  fprintf('Creating Generator Matrix ...\n');
  
  for i=1:k
    y = zeros(1,n);  
    y(k+1-i) = 1;

    y_gf = gf(y, m);
    % divide X^(n-k-i) by generator polynomial using function "deconv"
    % X^(n-k-i) => y_gf; generator polynomial => gf(gen_poly.x, m);
    [~,r] = deconv(y_gf, gf(gen_poly.x, m));
    % extract remainder
    remd(i,:) = r.x;
  end

  % Create Generator Matrix in systematic form
  G_matrix_full = cat(2, remd(:,k+1:n), eye(k));

  % Shortening procedure to align the number of source bits to "n_data_bits"
  G_matrix = G_matrix_full(1:k-n_delete, 1:n-n_delete);

  H_parity = G_matrix(:,1:n-k)';

  % Parity-check matrix n-k x n for BCH-DEC code
  H = cat(2, eye(n-k), H_parity);

end