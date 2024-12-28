% Author : Miroslav Marinkovic
% E-mail : mmarinkovic78d@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roots = roots_calculation(lambda, s_bch)

% roots calculation utilizes lookup tables (LUTs) to replace the Chien search, for t = 1, 2, 3;

% Reference literature:
% Xinmiao Zhang: "VLSI Architecture for Modern Error-Correcting Codes",
% page 42-49, 2016.
% Xinmiao Zhang, Michael O'Sullivan: "Ultra-Compressed Three-Error-Correcting BCH Decoder"
% IEEE ISCAS, 2018.

s = s_bch;
clear s_bch;
n = s.n;
t = s.t;
n_max = s.n_max;

% Create function handle
inv_gf = @(x) s.inv_x(uint32(x));

switch t
  
  case 1
    b = 0;
       
  case 2
    roots_two = @(x)s.roots_two(x,:);
    % Calculate variable a and c
    if (lambda(2) ~= 0)
      a = gf_mul(lambda(2), lambda(2), s);
      c = gf_mul(lambda(3), inv_gf(a), s);
    end
    if (lambda(3) ~= 0 && lambda(2) ~= 0)
      % Look-up table (LUT) for roots followed by inverse coordinate
      % transformation (to x)
      x = gf_mul(roots_two(c), lambda(2), s);
    end
    b = lambda(3);
    
  case 3
    % Create function handle
    roots_three = @(x)s.roots_three(x,:);
    roots_three_exc = @(x)s.roots_three_exception(x,:);
    square_root_gf = @(x) s.square_root_x(uint32(x));

    % Calculate variable a and b
    a = bitxor(gf_mul(lambda(2), lambda(2), s), lambda(3));
    b = bitxor(gf_mul(lambda(2), lambda(3), s), lambda(4));

    if (a ~= 0 && b ~= 0)
      square_a = gf_mul(a, a, s);
      cube_a = gf_mul(square_a, a, s);
      square_root_cube_a = square_root_gf(cube_a);
      c = gf_mul(b, inv_gf(square_root_cube_a), s);
      % Look-up table (LUT) for roots followed by inverse coordinate
      % transformation (to x)
      y = gf_mul(roots_three(c), square_root_gf(a), s);
      x = bitxor(y, lambda(2));
    elseif b ~= 0
      % Look-up table (LUT) for roots followed by inverse coordinate
      % transformation (to x)
      x = bitxor(uint32(roots_three_exc(b)), lambda(2));
    end
    
  otherwise
    error('error-correction capability is not supported');
end
    
if (b ~= 0) % double or triple error
  x = x(x > 0);
  roots = zeros(1, length(x));
  % if errors are correctable, roots belong to interval [1, n];
  % if errors are uncorrectable, roots can belong to interval [n+1, n_max];
  for i = 1:length(x)
    roots(i) = find(x(i) == s.alpha_power(1:n_max));
  end
else % single bit error, use lambda(2) to find a root
  roots = find(lambda(2) == s.alpha_power(1:n_max));
end

if isempty(roots)
  roots = 0;
else
  % if errors are correctable, roots belong to interval [1, n];
  % if errors are uncorrectable, roots can be negative if n < n_max
  roots = n + 1 - roots;
end

end

