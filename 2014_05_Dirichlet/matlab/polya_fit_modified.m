function a = polya_fit_modified(data, a)
% POLYA_FIT    Maximum-likelihood Dirichlet-multinomial (Polya) distribution.
%
% POLYA_FIT(data) returns the MLE (a) for the matrix DATA.
% Each row of DATA is a histogram of counts.
% POLYA_FIT(data,a) provides an initial guess A to speed up the search.
%
% The Polya distribution is parameterized as
%  p(x) = (Gamma(sum_k a_k)/prod_k Gamma(a_k)) prod_k Gamma(x_k+a_k)/Gamma(a_k)
%
% The algorithm is Newton iteration, described in
% "Estimating a Dirichlet distribution" by T. Minka.

% Written by Tom Minka

show_progress = 0;

if nargin < 2
  a = polya_moment_match(data);
end

ok = (col_sum(data) > 0);
if ~all(ok)
  a(ok) = polya_fit(data(:,ok), a(ok));
  return
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


% Newton-Raphson
old_e = sum(polya_logProb(a, data));
lambda = 0.1;
e = [];
for iter = 1:100
  if sum(a) == 0
    break
  end
  g = gradient2(a, U, V);
  abort = 0;
  % Newton iteration
  % loop until we get a nonsingular hessian matrix
  while(1)
    hg = 0 - hessian_times_gradient2(a, U, V, g, lambda);
    if all(hg < a)
      e(iter) = sum(polya_logProb(a-hg, data));
      if(e(iter) > old_e)
        old_e = e(iter);
        a = a - hg;
        lambda = lambda/10;
        break
      end
    end
    lambda = lambda*10;
    if lambda > 1e+6
      abort = 1;
      break
    end
  end
  if abort
    %disp('Search aborted')
    e(iter) = old_e;
    break
  end
  a(find(a < eps)) = eps;
  if max(abs(g)) < 1e-16
    break
  end
  if show_progress & rem(iter,5) == 0
    plot(e)
    drawnow
  end
end
if show_progress 
  disp(['gradient at exit = ' num2str(max(abs(g)))])
  plot(e)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = gradient2(a, U, V)

K = rows(U);
M = cols(U);

termToSubtract = 0;
for i = 1:M
  termToSubtract = termToSubtract + (V(i) / (sum(a) + i - 1));
end

g = zeros(1, K);
for k = 1:K
  for i = 1:M
    g(k) = g(k) + U(k, i) / (a(k) + i - 1);
  end
end

for k = 1:K
  g(k) = g(k) - termToSubtract;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hg = hessian_times_gradient2(a, U, V, g, lambda)

K = rows(U);
M = cols(U);

z = 0;
for i = 1:M
  z = z + V(i) / (sum(a) + i - 1)^2;
end

q = zeros(1, K);
for k = 1:K
  for i = 1:M
    q(k) = q(k) - U(k, i) / (a(k) + i - 1)^2;
  end
end

q = q - lambda;
q = 1./q;
b = sum(g .* q)/(1/z + sum(q));
hg = (b - g).*q;
