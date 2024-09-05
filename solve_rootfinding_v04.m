% Copyright contributors to the Hybrid-BPDN-solver project

function [x,r,g,data] = solve_rootfinding_v04(A,b,sigma,options)


t0 = tic();


if (isnumeric(A))
   [m,n] = size(A);
else
   info = A(0,0);
   m = info{1};
   n = info{2};
end

r = b; x = zeros(n,1); tau = 0;

if (isnumeric(A))
   g = A'*b;
else
   g = A(b,2);
end

data = struct();
data.info = {};
data.runtime = [];
data.iterations = 0;

iteration = 0; runtime = 0; innerIter = 0; innerIterH = 0;
info = struct('iter',0');
while (1)

   gNorm = norm(g,Inf);
   rNorm = sqrt(r'*r);

   time = toc(t0);
   data.runtime(end+1) = time;
   data.iterations(end+1) = info.iter;
   data.info{end+1} = info;

   fprintf('   %2d Tau = %f  rNorm = %f  sigma = %f  %.2s\n',iteration, tau,rNorm,sigma,time);
   t0 = tic();

   if (abs(rNorm - sigma) / max(1e-3,sigma) < 1e-5)
      break;
   end
   if (iteration == 20)
      break
   end

   % Update tau
   tau = tau + rNorm * (rNorm - sigma) / gNorm;

   [x,r,g,info] = solver_v04( A, b, [], 0, tau, [], x, options);
   iteration = iteration + 1;
end

data.tau = tau;
