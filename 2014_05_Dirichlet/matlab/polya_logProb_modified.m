function p = polya_logProb_modified(a, U, V)
% POLYA_LOGPROB   Dirichlet-multinomial (Polya) distribution.
%
% POLYA_LOGPROB(a,data) returns a vector containing the log-probability of 
% each histogram in DATA, under the Polya distribution with parameter A.
% U, V is a precomputed matrix and vector.
% If A is a row vector, then the histograms are the rows, otherwise columns.

K = rows(U);
M = cols(U);

if any(a < 0)
  p = -Inf;
  return
end
row = (rows(a) == 1);

p = 0;
for k = 1:K
  for i = 1:M
    p = p + U(k, i) * log(a(k) + i - 1);
  end
end

s = sum(a);
for i = 1:M
  p = p - V(i) * log(s + i - 1);
end
