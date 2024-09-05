% Copyright contributors to the Hybrid-BPDN-solver project

function [A,b,x0,sigma] = generate_problem_correlated(problem, m,n,k,mu,index)

s = RandStream('mt19937ar','Seed',3*m+17*n + k + 1872*index);
A = randn(s,m,n);
A = randomWalkMatrix(A,mu);

% Determine the support
v = randn(s,n,1);
[v,idx] = sort(v);

x0 = zeros(n,1);

switch (problem)
   case 1
      x0(idx(1:k)) = randn(s,k,1);      % Gaussian
      b = A*x0;
      sigma = [];

   case 2
      x0(idx(1:k)) = 2*rand(s,k,1) - 1; % Uniform
      b = A*x0;
      sigma = [];

   case 3
      x0(idx(1:k)) = sign(randn(s,k,1)); % Sign
      b = A*x0;
      sigma = [];

   case 4
      x0(idx(1:k)) = sign(randn(s,k,1)); % Sign + 1% noise
      b = A*x0;
      sigma = 0.01 * norm(b,2);
      z = randn(s,m,1);
      b = b + (sigma / norm(z,2)) * z;

   case 5
      x0(idx(1:k)) = sign(randn(s,k,1)); % Sign + 2% noise
      b = A*x0;
      sigma = 0.02 * norm(b,2);
      z = randn(s,m,1);
      b = b + (sigma / norm(z,2)) * z;

   case 6
      x0(idx(1:k)) = sign(randn(s,k,1)); % Sign + 5% noise
      b = A*x0;
      sigma = 0.05 * norm(b,2);
      z = randn(s,m,1);
      b = b + (sigma / norm(z,2)) * z;

   case 7
      x0(idx(1:k)) = sign(randn(s,k,1)); % Sign + 10% noise
      b = A*x0;
      sigma = 0.10 * norm(b,2);
      z = randn(s,m,1);
      b = b + (sigma / norm(z,2)) * z;

   otherwise
      error('Invalid problem index.');
end
