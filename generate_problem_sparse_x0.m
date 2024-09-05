% Copyright contributors to the Hybrid-BPDN-solver project

function [A,b,x0] = generate_problem_sparse_x0(m,n,k,type,idx)

s = RandStream('mt19937ar','Seed', 17*m + 131*n + type + idx*719);

A = randn(s, m,n) / sqrt(m);

% Determine the support
[v,idx] = sort(randn(s,1,n));
idx = idx(1:k);

% Initialize x0
x0 = zeros(n,1);

switch (type)
   case 1
       x0(idx) = sign(randn(s,1,k));

   case 2
       x0(idx) = 2 * rand(s,1,k) - 1;

   case 3
       x0(idx) = randn(s,1,k);

   otherwise
      error('Invalid problem type');
end

% Set the measurements
b = A*x0;
