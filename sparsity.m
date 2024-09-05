% Copyright contributors to the Hybrid-BPDN-solver project

function s = sparsity(x)

x = sort(abs(x),'descend');
s = zeros(1,9);
c = cumsum(x) / sum(x);
e = sqrt(cumsum(x.^2) / sum(x.^2));

s(1) = nnz(x);
s(2) = find(c > 0.980,1,'first');
s(3) = find(c > 0.990,1,'first');
s(4) = find(c > 0.999,1,'first');
s(5) = find(e > 0.980,1,'first');
s(6) = find(e > 0.990,1,'first');
s(7) = find(e > 0.999,1,'first');
s(8) = sum(x >= 1e-4*x(1));
s(9) = sum(x >= 1e-6*x(1));
