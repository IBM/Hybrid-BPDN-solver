% Copyright contributors to the Hybrid-BPDN-solver project

function [x,r,g,data] = solve_rootfinding_v05(A,b,sigma,options)
% Changes from v04: Maximum number of root finding iterations is now 1000

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
   if (iteration == 10000)
      break
   end

   % Update tau
   deltaTau = rNorm * (rNorm - sigma) / gNorm;
   if ((deltaTau <= 0) || (iteration == 0))
      multiplier = 1;
   else
      % See how many steps we can take while remaining feasible
      % This applies mostly to low-accuracy solves

      % Note that the last dual point may be (far) below fDualMax
      % so we use info.fDual instead of compute the gap based on the
      % current residual value r.
      tauIncr = options.optTol * max(info.f, options.optTolMinF);
      tauIncr = tauIncr - (info.f - info.fDual);
      tauIncr = tauIncr / gNorm;
      multiplier = max(1,ceil(tauIncr / deltaTau)); % Jump
      if (multiplier > 1)
         fprintf('   Jumping %d . . .\n', multiplier);
      end
   end
   tau = tau + multiplier * deltaTau;

   [x,r,g,info] = solver_v05( A, b, [], 0, tau, [], x, options);
   iteration = iteration + 1;
end

data.tau = tau;
