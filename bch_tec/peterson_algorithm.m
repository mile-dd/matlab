% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lambda = peterson_algorithm(syndrome, s_bch)

% Peterson痴 Algorithm to calculate error locator polynominal,
% for t = 1, 2, 3;

s = s_bch;
clear s_bch;

t = s.t; % error correction capability

% Create function handle
inv_gf = @(x) s.inv_x(uint32(x));

% Initialization of error locator polynominal
lambda = uint32(zeros(1, t+1));

switch t
  case 1
    lambda(1) = 1;
    lambda(2) = syndrome(1);
    
  case 2
    lambda(1) = 1;
    lambda(2) = syndrome(1);
    
    if syndrome(1) ~= 0      
      s1_square = gf_mul(syndrome(1), syndrome(1), s);
      lambda(3) = bitxor(s1_square, gf_mul(syndrome(3), inv_gf(syndrome(1)), s));
    end
    
  case 3
    s1_square = gf_mul(syndrome(1), syndrome(1), s);
    s1_cube = gf_mul(s1_square, syndrome(1), s);
    s1_cube_plus_s3 = bitxor(s1_cube, syndrome(3));
    
    lambda(1) = 1;
    lambda(2) = syndrome(1);
    
    if s1_cube ~= syndrome(3) % two or three errors
      s1_square_s3_plus_s5 = bitxor(gf_mul(s1_square, syndrome(3), s), syndrome(5));
      lambda(3) = gf_mul(s1_square_s3_plus_s5, inv_gf(s1_cube_plus_s3), s);
      lambda(4) = bitxor(s1_cube_plus_s3, gf_mul(syndrome(1), lambda(3), s));
    end
    
  otherwise
    error('error-correction capability t > 3 is not supported');
end

end

