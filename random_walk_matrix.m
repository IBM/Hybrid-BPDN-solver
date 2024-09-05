% Copyright contributors to the Hybrid-BPDN-solver project

function A = random_walk_matrix(A,mu) 

% Mu is the mutual coherence between successive cols
n = size(A,2);

% Find alpha
alphaMin = 0; alphaMax = 1;
for i=1:30
   alpha = (alphaMin + alphaMax) / 2;
   newMu = alpha / sqrt(alpha^2 + (1-alpha)^2);
   if (newMu > mu)
      alphaMax = alpha;
   else
      alphaMin = alpha;
   end
end

A = full(A * spdiags(1./sqrt(sum(A.^2,1)'),0,n,n));

for i=2:n
   % Two unit-norm columns
   a = A(:,i-1);
   b = A(:,i);
   v = b - (a'*b) * a;
   v = v / norm(v,2);

   c = alpha * a + (1-alpha) * v;
   c = c / norm(c,2);

   A(:,i) = c;
end
