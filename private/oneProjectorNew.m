function [x,lambda] = oneProjectorNew(b,d,tau)
% ONEPROJECTOR  Projects b onto the weighted one-norm ball of radius tau
%
%    [X,lambda] = ONEPROJECTOR(B,TAU) returns the orthogonal projection
%    of the vector b onto the one-norm ball of radius tau. The return
%    vector X which solves the problem
%
%            minimize  ||b-x||_2  st  ||x||_1 <= tau.
%               x
%
%    [X,lambda] = ONEPROJECTOR(B,D,TAU) returns the orthogonal
%    projection of the vector b onto the weighted one-norm ball of
%    radius tau, which solves the problem
%
%            minimize  ||b-x||_2  st  || Dx ||_1 <= tau.
%               x
%
%    If D is empty, all weights are set to one, i.e., D = I.
%
%    In both cases, the return value lambda gives the soft-thresholding
%    value.
%
% See also spgl1.


% Check arguments
if nargin < 2
  error('The oneProjector function requires at least two parameters');
end
if nargin < 3
  tau = d;
  d   = [];
end

% Check weight vector
if isempty(d), d = 1; end;

if ~isscalar(d) && ( length(b) ~= length(d) )
  error('Vectors b and d must have the same length');
end

% Quick return for the easy cases.
if ((isscalar(d)) && (d == 0))
   x     = b;
   lamda = 0;
   return
end

% Get sign of b and set to absolute values
s = sign(b);
b = abs(b);

% Perform the projection
if isscalar(d)
   w = b;
   wNnz  = nnz(w);
   wNorm = sum(w);
   wMax  = max(w);
   
   for i=1:1000 % Large enough number to work in vast majority of cases
      wNormPrev = wNorm;
   
      lambda = (wNorm-tau) / wNnz;
      if (lambda > (1-1e-8) * wMax)
         lambda = (1-1e-8) * wMax;
      end
   
      wNew  = max(0,w - lambda);
      wNorm = norm(wNew,1);

      if (wNorm < tau)
         w = max(0,w - lambda*0.5);
         wNorm = norm(w,1);
      else
         w = wNew;
      end

      wNnz = nnz(w);      
      wMax = max(w);
   
      if ((wNorm <= (1+eps)*tau) || (wNorm == wNormPrev))
         break;
      end
   end
   
   x = w;
else
   %d   = abs(d);
   %idx = find(d > eps); % Get index of all non-zero entries of d
   %x   = b;             % Ensure x_i = b_i for all i not in index set idx
   %[x(idx),lambda] = oneProjectorMex(b(idx),d(idx),tau);
   error('TODO');
end

% Restore signs in x
x = x.*s;





