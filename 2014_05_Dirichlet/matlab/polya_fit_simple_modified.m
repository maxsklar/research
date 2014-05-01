function a = polya_fit_simple(data,a)
% POLYA_FIT_SIMPLE   Maximum-likelihood Polya distribution.
%
% Same as POLYA_FIT but uses the simple fixed-point iteration described in
% "Estimating a Dirichlet distribution" by T. Minka. 

show_progress = 0;

if nargin < 2
  a = polya_moment_match(data);
end
sdata = row_sum(data);

% Finding the U-matrix and v-vector
M = max(sdata);
K = size(data, 2);
U = zeros(K, M);
V = zeros(1, M);

for n = 1:rows(data)
  for k = 1:K
    for count = 1:data(n, k)
      U(k, count) = U(k, count) + 1;
    end
  end
  for count = 1:sdata(n)
    V(count) = V(count) + 1;
  end
end

% fixed-point iteration
[N,K] = size(data);
for iter = 1:1000
  old_a = a;
  sa = sum(a);
  g = numerator(a, U);
  h = denominator(sa, V);
  
  a = a .* g ./ h;
  if show_progress
    e(iter) = sum(polya_logProb(a, data));
  end
  if max(abs(a - old_a)) < 1e-6
    break
  end
  if show_progress & rem(iter,10) == 0
    plot(e)
    drawnow
  end
end  
if show_progress
  plot(e)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function n = numerator(a, U)

K = rows(U);
M = cols(U);

n = zeros(1, K);
for k = 1:K
  for i = 1:M
    n(k) = n(k) + U(k, i) / (a(k) + i - 1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = denominator(sa, V)

M = cols(V);

d = 0;
for i = 1:M
  d = d + V(i) / (sa + i - 1);
end
