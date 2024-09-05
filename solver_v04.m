% Copyright contributors to the Hybrid-BPDN-solver project

function [x,r,g,info] = solver_v04( A, b, c, mu, tau, sigma, x, options )

% CHANGES TO:
% v03 USE A FAST IMPLEMENTATION OF TRAJECTORY SEARCH
% v04 USE OPTTOLMINF AS INTENDED

%SPGL1  Solve basis pursuit, basis pursuit denoise, and LASSO
%
% [x, r, g, info] = solver_v01(A, b, c, mu, tau, sigma, x0, options)
%
% ---------------------------------------------------------------------
% Solve the augmented basis pursuit denoise (BPDN) problem
%
% (BPDN)   minimize ||x||_1
%
%          subject to  1/2*||Ax-b||_2^2 + c^Tx + (mu/2)*||x||_2^2 <= sigma^2,
%
% or the l1-constrained regularized least-squares problem
%
% (LASSO)  minimize  1/2*||Ax-b||_2^2 + c^Tx + (mu/2)*||x||_2^2
%
%          subject to  ||x||_1 <= tau.
% ---------------------------------------------------------------------
%
% INPUTS
% ======
% A        is an m-by-n matrix, explicit or an operator.
%          If A is a function, then it must have the signature
%
%          y = A(x,mode)   if mode == 1 then y = A x  (y is m-by-1);
%                          if mode == 2 then y = A'x  (y is n-by-1).
%
% b        is an m-vector.
% c        is an n-vector.
% mu       is a scalar >= 0.
% tau      is a nonnegative scalar; see (LASSO).
% sigma    if sigma != inf or != [], then spgl1 will launch into a
%          root-finding mode to find the tau above that solves (BPDN).
%          In this case, it's STRONGLY recommended that tau = 0.
% x0       is an n-vector estimate of the solution (possibly all
%          zeros). If x0 = [], then SPGL1 determines the length n via
%          n = length( A'b ) and sets  x0 = zeros(n,1).
% options  is a structure of options from spgSetParms. Any unset options
%          are set to their default value; set options=[] to use all
%          default values.
%
% OUTPUTS
% =======
% x        is a solution of the problem
% r        is the residual, r = b - Ax
% g        is the gradient, g = -A'r
% info     is a structure with the following information:
%          .tau     final value of tau (see sigma above)
%          .rNorm   two-norm of the optimal residual
%          .rGap    relative duality gap (an optimality measure)
%          .gNorm   Lagrange multiplier of (LASSO)
%          .stat    = 1 found a BPDN solution
%                   = 2 found a BP sol'n; exit based on small gradient
%                   = 3 found a BP sol'n; exit based on small residual
%                   = 4 found a LASSO solution
%                   = 5 error: too many iterations
%                   = 6 error: linesearch failed
%                   = 7 error: found suboptimal BP solution
%                   = 8 error: too many matrix-vector products
%          .time    total solution time (seconds)
%          .nProdA  number of multiplications with A
%          .nProdAt number of multiplications with A'
%
% OPTIONS
% =======
% Use the options structure to control various aspects of the algorithm:
%
% options.fid         File ID to direct log output
%        .verbosity   0=quiet, 1=some output, 2=more output.
%        .iterations  Max. number of iterations (default if 10*m).
%        .bpTol       Tolerance for identifying a basis pursuit solution.
%        .optTol      Optimality tolerance (default is 1e-4).
%        .decTol      Larger decTol means more frequent Newton updates.
%
% EXAMPLE
% =======
%   m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
%   p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
%   A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
%   b  = A*x0 + 0.005 * randn(m,1);
%   opts = spgSetParms('optTol',1e-4);
%   [x,r,g,info] = spgl1(A, b, 0, 1e-3, [], opts); % Find BP sol'n.
%
% AUTHORS
% =======
%  Ewout van den Berg (ewout78@cs.ubc.ca)
%  Michael P. Friedlander (mpf@cs.ubc.ca)
%    Scientific Computing Laboratory (SCL)
%    University of British Columbia, Canada.
%
% BUGS
% ====
% Please send bug reports or comments to
%            Michael P. Friedlander (mpf@cs.ubc.ca)
%            Ewout van den Berg (ewout78@cs.ubc.ca)

% 15 Apr 07: First version derived from spg.m.
%            Michael P. Friedlander (mpf@cs.ubc.ca).
%            Ewout van den Berg (ewout78@cs.ubc.ca).
% 17 Apr 07: Added root-finding code.
% 18 Apr 07: sigma was being compared to 1/2 r'r, rather than
%            norm(r), as advertised.  Now immediately change sigma to
%            (1/2)sigma^2, and changed log output accordingly.
% 24 Apr 07: Added quadratic root-finding code as an option.
% 24 Apr 07: Exit conditions need to guard against small ||r||
%            (ie, a BP solution).  Added test1,test2,test3 below.
% 15 May 07: Trigger to update tau is now based on relative difference
%            in objective between consecutive iterations.
% 15 Jul 07: Added code to allow a limited number of line-search
%            errors.
% 23 Feb 08: Fixed bug in one-norm projection using weights. Thanks
%            to Xiangrui Meng for reporting this bug.
% 26 May 08: The simple call spgl1(A,b) now solves (BPDN) with sigma=0.

%   spgl1.m
%   $Id: spgl1.m 1407 2009-06-30 20:00:54Z ewout78 $
%
%   ----------------------------------------------------------------------
%   This file is part of SPGL1 (Spectral Projected-Gradient for L1).
%
%   Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
%   Department of Computer Science, University of British Columbia, Canada.
%   All rights reserved. E-mail: <{ewout78,mpf}@cs.ubc.ca>.
%
%   SPGL1 is free software; you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as
%   published by the Free Software Foundation; either version 2.1 of the
%   License, or (at your option) any later version.
%
%   SPGL1 is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
%   Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with SPGL1; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
%   USA
%   ----------------------------------------------------------------------
REVISION = '$Revision: 1017 $';
DATE     = '$Date: 2008-06-16 22:43:07 -0700 (Mon, 16 Jun 2008) $';
REVISION = REVISION(11:end-1);
DATE     = DATE(35:50);

% VERSION 3: Use Mex implementation of Product B.
% VERSION 4: Fixes bug in version 3: we need to check if the projected
%            gradient lies in the tangent wedge of the face before we
%            can apply the quasi-Newton step on the face. Failure to do
%            this can lead to a minimum in the relative interior of the
%            face that is not a global minimum.
% VERSION 5: - Fixed bugs in self-projection cone
%            - Fixed bug in faces, we need to look at support and one-norm

tic; % Start your watches!
m = length(b);

%----------------------------------------------------------------------
% Check arguments.
%----------------------------------------------------------------------
if ~exist('options','var'), options = []; end
if ~exist('x',      'var'), x       = []; end
if ~exist('c',      'var'), c       = []; end
if ~exist('mu',     'var'), mu      =  0; end
if ~exist('tau',    'var'), tau     = []; end
if ~exist('sigma',  'var'), sigma   = []; end

if nargin < 2 || isempty(b) || isempty(A)
   error('At least two arguments are required');
elseif isempty(tau) && isempty(sigma)
   tau = 0;
   sigma = 0;
   singleTau = false;
elseif isempty(sigma) % && ~isempty(tau)  <-- implied
   singleTau = true;
else
   if isempty(tau)
      tau = 0;
   end
   singleTau = false;
end

if (all(c==0))
   c = [];
end

if (tau   < 0), error('Tau must be greater than or equal to 0.');   end
if (sigma < 0), error('Sigma must be greater than or equal to 0.'); end
if (mu    < 0), error('Mu must be greater than or equal to 0.');    end

if ((~singleTau) && (~isempty(c) || mu > 0))
   error('Basis-pursuit formulation is currently supported only for c=0 and mu=0.');
end

%----------------------------------------------------------------------
% Grab input options and set defaults where needed.
%----------------------------------------------------------------------
fid          = option('fid',            1);    % File ID for output
logLevel     = option('verbosity',      2);    % Verbosity level
history      = option('history',    false); % Set flag to record iteration information

% Stopping criteria
bpTol        = option('bpTol',      1e-06); % Tolerance for basis pursuit solution
lsTol        = option('lsTol',      1e-06); % Least-squares optimality tolerance
optTol       = option('optTol',     1e-04); % Optimality tolerance
optTolRel    = option('optTolRel',  false); %   False uses min(optTolMinF,f)
optTolMinF   = option('optTolMinF',     1); %   True always uses f (same as optTolMinF = -Inf)
decTol       = option('decTol',     1e-04); % Req'd rel. change in primal obj. for Newton
maxIts       = option('iterations',  10*m); % Max number of iterations
maxMatvec    = option('maxMatVec',    Inf); % Maximum matrix-vector multiplies allowed

% Line search
stepMin      = option('stepMin',    1e-16); % Minimum spectral step
stepMax      = option('stepMax',    1e+05); % Maximum spectral step
nPrevVals    = option('nPrevVals',      3); % Number previous func values for linesearch
lsType       = option('lsType',         1); % Type of line search

% Quasi-Newton
qnType       = option('qnType',         0); % Quasi-Newton type
qnAllowFull  = option('qnAllowFull',    0); % Allows quasi-Newton on full dimension
qnStrict     = option('qnStrict',       1); % Allows quasi-Newton on dimension <= m only
qnWolfe      = option('qnWolfe',    false); % Enforce Wolfe line-search condition (beta = 0.99)
qnResetHist  = option('qnResetHist',false); % Reset function value history after step

% Dual type
dualType     = option('dualType',       0); % 0 = original, 1 = new dual
dualBest     = option('dualBest',    true); % true = use best, false = use latest


% Root finding
rfType       = option('rfType',         1); % Root-finding type

maxLineErrors = 10;     % Maximum number of line-search failures.

% Exit conditions (constants).
EXIT_ROOT_FOUND    = 1;
EXIT_BPSOL_FOUND   = 2;
EXIT_LEAST_SQUARES = 3;
EXIT_OPTIMAL       = 4;
EXIT_ITERATIONS    = 5;
EXIT_LINE_ERROR    = 6;
EXIT_SUBOPTIMAL_BP = 7;
EXIT_MATVEC_LIMIT  = 8;

% Line search types
LS_CLASSIC           = 1; % Curvilinear and backtracking projected
LS_CURVE_EXACT_PROJ  = 2; % Curvilinear and exact from projected
LS_EXACT_PROJ        = 3; % Exact from projected only
LS_BACKTRACK_PROJ    = 4; % Backtracking from projected only
LS_TRAJECTORY        = 5; % First local minimum on trajectory
LS_TRAJECTORY_GLOBAL = 6; % Global minimum on trajectory
LS_CLASSIC_PLUS      = 7; % Same as classic, but with exact from projected
                          % added as a final back-up option

% Root-finding types
RF_SPGL1            = 1;  % Classic mode
RF_STRICT           = 2;

% Quasi-Newton types
QN_NONE             = 0;  % Disable quasi-Newton search-direction
QN_CONSECUTIVE      = 1;  % Enable when support stays the same
QN_PARENT_FAST      = 2;  % Enable when support same or grows
QN_PARENT_SLOW      = 3;  % Same as previous, use one iteration later

if (option('classic',false) == true)
   lsType  = LS_CLASSIC;
   lsFirst = false;
   qnType  = QN_NONE;
   rfType  = RF_SPGL1;
   classicMode = true;
else
   classicMode = false;
end

if (qnType ~= QN_NONE)
   lbfgsHist = 8;
end


%----------------------------------------------------------------------
% Initialize local variables.
%----------------------------------------------------------------------
iter          = 0;
nProdA        = 0;
nProdAt       = 0;
lastFv        = -inf(nPrevVals,1);  % Last m function values.
nLineTot      = 0;                  % Total no. of linesearch steps.
printTau      = false;
nNewton       = 0;
bNorm         = norm(b,2);
stat          = false;
timeProject   = 0;
timeMatProd   = 0;
stepG         = 1;                  % Step length for projected gradient.
testUpdateTau = 0;                  % Previous step did not update tau

if (history)
   historyTime = zeros(1,maxIts+1);
   historyFun  = zeros(1,maxIts+1);
   historyDual = zeros(1,maxIts+1);
   historyGap  = zeros(1,maxIts+1);
   historyQN   = false(1,maxIts+1);
   historyAprod= zeros(2,maxIts+1);
end

% Determine initial x and vector length n
explicit = ~(isa(A,'function_handle'));
if isempty(x)
   if isnumeric(A)
      n = size(A,2);
   else
      x = Aprod(b,2);
      n = length(x);
   end
   x = zeros(n,1);
else
   n = length(x);
end

% TODO: weights
weights     = 1;
weightsFull = weights * ones(n,1);


% Check for complex numbers
if ((qnType ~= QN_NONE) && (~isreal(x) || ~isreal(b) || ~isreal(c)))
   warning('Disabled quasi-Newton for complex problem.');
   qnType = QN_NONE; % Disable quasi-Newton for complex problems
end


%----------------------------------------------------------------------
% Log header.
%----------------------------------------------------------------------
printf('\n');
printf(' %s\n',repmat('=',1,80));
printf(' SPGL1  v.%s (%s)\n', REVISION, DATE);
printf(' %s\n',repmat('=',1,80));
printf(' %-22s: %8i %4s'   ,'No. rows'          ,m       ,'');
printf(' %-22s: %8i\n'     ,'No. columns'       ,n          );
printf(' %-22s: %8.2e %4s' ,'Initial tau'       ,tau     ,'');
printf(' %-22s: %8.2e\n'   ,'Two-norm of b'     ,bNorm      );
printf(' %-22s: %8.2e %4s' ,'Optimality tol'    ,optTol  ,'');
if singleTau
   printf(' %-22s: %8.2e\n'  ,'Target one-norm of x'  ,tau  );
else
   printf(' %-22s: %8.2e\n','Target objective'  ,sigma      );
end
printf(' %-22s: %8.2e %4s' ,'Basis pursuit tol' ,bpTol   ,'');
printf(' %-22s: %8i\n'     ,'Maximum iterations',maxIts     );
printf('\n');
if singleTau
   logB = ' %5i  %13.7e  %13.7e  %9.2e  %6.1f  %7.1e/%7.1e  %6i  %s';
   logH = ' %5s  %13s  %13s  %9s  %6s  %6s  %7s/%7s   %6s\n';
else
   logB = ' %5i  %13.7e  %13.7e  %9.2e  %6.1f  %7.1e/%7.1e  %6i  %s';
   logH = ' %5s  %13s  %13s  %9s  %6s  %6s  %7s/%7s   %6s\n';
end


%----------------------------------------------------------------------
% Get started.
%----------------------------------------------------------------------

% Quick exit if sigma >= ||b||.  Set tau = 0 to short-circuit the loop.
if (~isempty(sigma) && (bNorm <= sigma) && (isempty(c)))
   printf('W: sigma >= ||b||.  Exact solution is x = 0.\n');
   tau = 0;  singleTau = true;
end

% Project the starting point and evaluate function and gradient.
x = project(x,tau,weights);
r = b - Aprod(x,1);  % r = b - Ax
g =   - Aprod(r,2);  % g = -A'r
f = (r'*r) / 2;

if (~isempty(c))
   g = g + c;
   f = f + c'*x;
end
if (mu > 0)
   g = g + mu * x;
   f = f + (mu/2) * (x'*x);
end

% Set current maximum dual objective
fDualMax = -Inf;

% Required for nonmonotone strategy.
lastFv(1) = f;
fBest     = f;
xBest     = x;
fOld      = f;
xOld      = x;

% Compute projected gradient direction and initial steplength.
dx     = project(x - g, tau, weights) - x;
dxNorm = norm(dx,inf);
if dxNorm < (1 / stepMax)
   bbScaling = stepMax;
else
   bbScaling = min( stepMax, max(stepMin, 1/dxNorm) );
end

% Maintain support
if (qnType ~= QN_NONE)
   xAbs     = abs(x);
   xNorm1   = sum(xAbs);
   if ((abs(xNorm1 - tau) / max(1,tau)) < 1e-8)
      support = false(n,1); % Interior
   else
      support  = (xAbs >= 1e-9); % Boundary
   end
   supportChange = 0;
else
   support = [];
   supportChange = [];
end

iterH = 0;    % Number of quasi-Newton iterations
H     = [];   % Current estimate of the Hessian
sqrt1 = [];   % Intermediate square-root values #1
sqrt2 = [];   % Intermediate square-root values #2
flagUseHessian = false;

statusMsg = '';

%----------------------------------------------------------------------
% MAIN LOOP.
%----------------------------------------------------------------------
while 1
    % =================================================================
    % Test exit conditions.
    % =================================================================

    % Compute dual objectives
    gNorm = NormL1_dual(-g,weights); % ||A^Tr - mu*x - c||

    % First way of computing the duality gap
    if ((mu == 0) || (dualType == 0))
       fDual = r'*b - tau*gNorm - (r'*r)/2;
       if (mu > 0)
          fDual = fDual - (mu/2)*(x'*x);
       end
    end

    % Second way of computing the duality gap when mu > 0
    if ((mu > 0) && (dualType == 1))
       % Determine z = |A^Tr - c^Tx| = |-g + mu*x|
       z = abs(mu * x - g);

       % Solve: Minimize  tau * lamba + (1/2*mu)||[z - \lambda*w]_+||_2^2
       %        lambda>=0
       [lambdaStar, objValue] = find_lambda_star(z,weightsFull,tau,mu);

       % Compute the duality objective
       fDual = r'*b - (r'*r) / 2 - objValue;
    end

    if (dualBest)
       if (fDual > fDualMax)
          fDualMax = fDual;
       else
          fDual = fDualMax;
       end
    end
    gap = f - fDual;

    % Compute quantities used in stopping criterion
    rNorm   = norm(r, 2);

    if (optTolRel)
       rGap = abs(gap) / f;
    else
       rGap = abs(gap) / max(optTolMinF,f);
    end

    aError1 = rNorm - sigma;
    aError2 = f - sigma^2 / 2;
    rError1 = abs(aError1) / max(1,rNorm);
    rError2 = abs(aError2) / max(1,f);

    % Record history information
    if (history)
       historyTime(iter+1)    = toc;
       historyFun(iter+1)     = f;
       historyDual(iter+1)    = fDual;
       historyAprod(1,iter+1) = nProdA;
       historyAprod(2,iter+1) = nProdAt;

       historyGap(iter+1) = gap;
    end

    % Single tau: Check if we're optimal.
    if singleTau
       if rGap <= optTol
          stat  = EXIT_OPTIMAL;
       end

    % Multiple tau: Check if found root and/or if tau needs updating.
    else
       if (rfType == RF_STRICT)
          testUpdateTau = false;
          stat = 0;

          if (rNorm < bpTol * bNorm)
             stat = EXIT_BPSOL_FOUND;
          end

          if (rGap <= decTol)
             if (rError1 < optTol)
                stat = EXIT_ROOT_FOUND;
             else
                testUpdateTau = true;
             end
          end
       end


       if (rfType == RF_SPGL1)
          % Test if a least-squares solution has been found
          if gNorm <= lsTol * rNorm
             stat = EXIT_LEAST_SQUARES;
          end

          if rGap <= max(optTol, rError2) || rError1 <= optTol
             % The problem is nearly optimal for the current tau.
             % Check optimality of the current root.
             test1 = rNorm       <=   bpTol * bNorm;
          %  test2 = gNorm       <=   bpTol * rNorm;
             test3 = rError1     <=  optTol;
             test4 = rNorm       <=  sigma;

             if test4, stat=EXIT_SUBOPTIMAL_BP;end  % Found suboptimal BP sol.
             if test3, stat=EXIT_ROOT_FOUND;   end  % Found approx root.
             if test1, stat=EXIT_BPSOL_FOUND;  end  % Resid minim'zd -> BP sol.
           % 30 Jun 09: Large tau could mean large rGap even near LS sol.
           %            Move LS check out of this if statement.
           % if test2, stat=EXIT_LEAST_SQUARES; end % Gradient zero -> BP sol.
          end

          testRelChange1 = (abs(f - fOld) <= decTol * f);
          testRelChange2 = (abs(f - fOld) <= 1e-1 * f * (abs(rNorm - sigma)));
          testUpdateTau  = ((testRelChange1 && rNorm >  2 * sigma) || ...
                            (testRelChange2 && rNorm <= 2 * sigma)) && ...
                            ~stat && ~testUpdateTau;

          if (~classicMode)
             if (rGap <= optTol)
                % The original version of SPGL1 can have a zero rGap
                % and still not update tau!
                testUpdateTau = true;
             end
          end
       end


       if (testUpdateTau && ~stat)
          % Update tau
          tauOld   = tau;
          tau      = max(0,tau + (rNorm * aError1) / gNorm);
          nNewton  = nNewton + 1;
          printTau = (abs(tauOld - tau) >= 1e-6 * tau); % For log only.
          if (tau < tauOld)
             % The one-norm ball has decreased. Project the current iterate
             % to ensure feasibility.
             x = project(x,tau,weights);

             % We have changed x and must recompute the residual, gradient
             % and objective value.
             r = b - Aprod(x,1);  % r = b - Ax
             g =   - Aprod(r,2);  % g = -A'r
             f = r'*r / 2;

             % The function history must also be reset otherwise there may
             % be no hope of satisfying the line search condition.
             lastFv = -inf(nPrevVals,1);  % Last m function values.
             lastFv(1) = f;

             % Update support
             if (qnType ~= QN_NONE)
                support  = (xAbs >= 1e-9); % Boundary
             end
          else
             if ((tau > tauOld) && (qnType ~= QN_NONE))
                support = false(n,1); % Interior
             end
          end

          % Reset Hessian
          H = []; flagUseHessian = false;

          % Reset status
          stat = 0;

          % Reset the current maximum dual objective
          fDualMax  = -Inf;
       end
    end

    % Too many its and not converged.
    if ~stat  &&  iter >= maxIts
        stat = EXIT_ITERATIONS;
    end


    % =================================================================
    % Print log, update history and act on exit conditions.
    % =================================================================
    if logLevel >= 2 || singleTau || printTau || iter == 0 || stat
       tauFlag = '              ';
       if printTau, tauFlag = sprintf(' %13.7e',tau);   end

       if singleTau
          if (mod(iter,50) == 0)
             printf(logH,'Iter','Objective','Relative Gap','gNorm','stepG','BB','||dx||','nnzX');
             printf('\n');
          end
          printf(logB,iter,f,rGap,gNorm,log10(stepG),bbScaling,norm(x-xOld),sum(support),statusMsg);
       else
          if (mod(iter,50) == 0)
             printf(logH,'Iter','Objective','Relative Gap','gNorm','stepG','BB','||dx||','nnzX');
             printf('\n');
          end
          printf(logB,iter,f,rGap,gNorm,log10(stepG),bbScaling,norm(x-xOld),sum(support),statusMsg);
          if printTau
             printf(' %s',tauFlag);
          end
       end

       printf('\n');
    end
    printTau = false;

    if stat, break; end % Act on exit conditions.


    %==================================================================
    %==================================================================
    %
    %           I T E R A T I O N S   B E G I N   H E R E
    %
    %==================================================================
    %==================================================================

    iter = iter + 1;
    xOld = x;  fOld = f;  gOld = g;  rOld = r; supportOld = support;

    statusMsg = '';

    try

       % ===================================
       % Try a quasi-Newton direction first
       % ===================================
       if (flagUseHessian)
          statusMsg = 'H'; iterH = iterH + 1;

          % --------------------------------------------------
          % Step 1. Get search direction
          % --------------------------------------------------
          if (~any(support))
             d = lbfgshprod(H,-g);
          else
             d = -g;

             % Project gradient onto the coefficient space
             dTrans = productBMex(signs .* d(support), 1, sqrt1, sqrt2);

             % If ||dTrans|| is tiny it is (near) orthogonal to the face
             if (norm(dTrans,2) <= 1e-10 * max(1,norm(d,2)))
                stat = EXIT_OPTIMAL;
             end

             % Get quasi-Newton search direction
             dQuasi = lbfgshprod(H,dTrans);

             % Convert back to global domain to get new search direction
             d = zeros(n,1);
             dSupport =  signs .* productBMex(dQuasi, 0, sqrt1, sqrt2);
             d(support) = dSupport;
          end

          % --------------------------------------------------
          % Step 2. Determine first non-zero entry to hit zero
          % --------------------------------------------------
          s1    = (((d < 0) & (x > 0)) | ((d > 0) & (x < 0)));
          if (any(s1))
             gamma = -max(x(s1) ./ d(s1));
          else
             gamma = +Inf;
          end

          w = Aprod(d,1);

          lnErr = false;

          % --------------------------------------------------
          % Step 3. Compute the optimal step length beta
          % --------------------------------------------------
          enumerator  = (w'*r);
          denominator = (w'*w);
          if (~isempty(c))
             enumerator = enumerator - (c'*d);
          end
          if (mu > 0)
             enumerator  = enumerator - mu * (x'*d);
             denominator = denominator + mu * (d'*d);
          end
          beta = enumerator / denominator;
          if (beta <= 1e-11)
             lnErr = true;
          else
             beta = min(beta,gamma);
          end

          % TODO: ADD SECOND CRITERION FOR WOLFE [Sufficient descent]
          %{
          if (qnWolfe)
             gtd = -((dTrans)'*dQuasi); % Same as g'*d
             if (gtd + beta * (w'*w) > 0.95*gtd) % 0.99
                lnErr = false;
             else
                lnErr = true;
             end
          else
             lnErr = false;
          end
          %}

          % --------------------------------------------------
          % Check if a suitable step length was found
          % --------------------------------------------------
          if ~lnErr
             x = xOld + beta * d;
             r = r - beta * w; % Avoid evaluating A*x
             f = (r'*r) / 2;

             stepG = beta;
             nLineTot = nLineTot + 1;

             % Reset function history
             if (qnResetHist)
                lastFv(:) = f;
             end

             % Record successful QN step
             if (history)
                historyQN(iter+1) = true;
             end
          end

       else
          % Do not use Hessian
          lnErr = true;
       end


       %---------------------------------------------------------------
       % Projected gradient step and linesearch.
       %---------------------------------------------------------------

       % -----------------------------------------------------------------
       % Line search #1: First local minimum along projection trajectory
       % -----------------------------------------------------------------
       if ((lnErr) && ...
           ((lsType == LS_TRAJECTORY) || (lsType == LS_TRAJECTORY_GLOBAL)))

          % Trajectory line search
          x = oneCauchyPointMex(A,b,x,-g,weightsFull,tau, lsType == LS_TRAJECTORY);
          nProdA = nProdA + 3;

          if (~isempty(x))
             r = b - Aprod(x,1);

             % Verify sufficient descent condition
             f = (r'*r) / 2;
             if (~isempty(c))
                f = f + c'*x;
             end
             if (mu > 0)
                f = f + (mu/2) * (x'*x);
             end

             if (f >= fOld - 1e-4 * (g' * (xOld - x)))
                % Sufficient descent not satisfied
                x = xOld;
                r = rOld;
                f = fOld;
                lnErr = true;
             else
                lnErr = false;
             end
          else
             x = xOld;
             r = rOld;
             f = fOld;
             lnErr = true;
          end

          statusMsg = [statusMsg,'T'];
       end


       % -----------------------------------------------------------------
       % Line search #2: Projected backtracking in which each trial point
       %                 is projected onto the one-norm ball
       % -----------------------------------------------------------------
       if ((lnErr) && ...
           ((lsType == LS_CLASSIC         ) || ...
            (lsType == LS_CLASSIC_PLUS    ) || ...
            (lsType == LS_CURVE_EXACT_PROJ) || ...
            (lsType == LS_TRAJECTORY      ) || ...
            (lsType == LS_TRAJECTORY_GLOBAL)))
          [f,x,r,nLine,stepG,lnErr] = ...
              spgLineCurvy(x,bbScaling*g,max(lastFv),@Aprod,b,c,mu,@project,tau,weights);
          nLineTot = nLineTot + nLine;

          if (lnErr)
             x = xOld;
             r = rOld;
             f = fOld;
          end

          statusMsg = [statusMsg,'C'];
       end

       % -----------------------------------------------------------------
       % Line search #3: Backtracking line search over segment between
       %                 projected point and the current x
       % -----------------------------------------------------------------
       if ((lnErr) && ...
           ((lsType == LS_CLASSIC          ) || ...
            (lsType == LS_CLASSIC_PLUS    ) || ...
            (lsType == LS_BACKTRACK_PROJ   ) || ...
            (lsType == LS_TRAJECTORY       ) || ...
            (lsType == LS_TRAJECTORY_GLOBAL)))

          %  Projected backtrack failed. Retry with feasible dir'n linesearch.
          dx   = project(x - bbScaling*g, tau, weights) - x;
          gtd  = g'*dx;
          [f,x,r,nLine,lnErr] = spgLine(f,x,dx,gtd,max(lastFv),@Aprod,b,c,mu);
          nLineTot = nLineTot + nLine;

          if (lnErr)
             x = xOld;
             r = rOld;
             f = fOld;
          end

          statusMsg = [statusMsg,'B'];
       end


       % -----------------------------------------------------------------
       % Line search #4: Exact minimum on segment between projected point
       %                 and the current x
       % -----------------------------------------------------------------
       if ((lnErr) && ...
           ((lsType == LS_CLASSIC_PLUS     ) || ...
            (lsType == LS_CURVE_EXACT_PROJ ) || ...
            (lsType == LS_EXACT_PROJ       ) || ...
            (lsType == LS_TRAJECTORY       ) || ...
            (lsType == LS_TRAJECTORY_GLOBAL)))
          %  Projected backtrack failed. Retry with feasible dir'n linesearch.
          p = project(x - bbScaling*g, tau, weights);
          v = p - x;
          w = Aprod(v,1);

          % Compute optimal beta -- Ignore any sufficient descent
          enumerator  = (w'*r);
          denominator = (w'*w);
          if (~isempty(c))
             enumerator = enumerator - (c'*v);
          end
          if (mu > 0)
             enumerator  = enumerator - mu * (x'*v);
             denominator = denominator + mu * (v'*v);
          end
          beta = enumerator / denominator;

          % Sufficient descent holds whenever gamma1 < 0.5

          % Check beta
          if (beta <= 1e-12)
             beta  = max(beta,0);
             lnErr = true;
          else
             beta  = min(beta,1);
             lnErr = false;
          end

          % Update iterate and objective value
          if (~lnErr)
             x = xOld + beta * v;
             Ax= Aprod(x,1);
             r = b - Ax;
             f = (r'*r) / 2;
             if (~isempty(c))
                f = f + c'*x;
             end
             if (mu > 0)
                f = f + (mu/2) * (x'*x);
             end

             stepG    = beta;
             nLineTot = nLineTot + 1;
          end

          statusMsg = [statusMsg,'P'];
       end

       % -----------------------------------------------------------------
       % Line search: Unable to find a suitable step size; scale down the
       %              Barzilai-Borwein scaling parameter
       % -----------------------------------------------------------------
       if ((lnErr) && ...
           ((lsType == LS_CLASSIC         ) || ...
            (lsType == LS_CLASSIC_PLUS    ) || ...
            (lsType == LS_EXACT_PROJ      ) || ...
            (lsType == LS_CURVE_EXACT_PROJ)))
          %  Failed again.  Revert to previous iterates and damp max BB step.
          if maxLineErrors <= 0
             stat = EXIT_LINE_ERROR;
          else
             stepMax = stepMax / 10;
             printf(['W: Linesearch failed with error %i. '...
                     'Damping max BB scaling to %6.1e.\n'],lnErr,stepMax);
             maxLineErrors = maxLineErrors - 1;
          end
       else
          if (lnErr), stat = EXIT_LINE_ERROR; end
       end

       absx = abs(x); % Used for support determination later
       xNorm1 = sum(absx); % norm(x,1);
       ensure(xNorm1 <= tau+optTol);


       %---------------------------------------------------------------
       % Update gradient and compute new Barzilai-Borwein scaling.
       %---------------------------------------------------------------
       if (~lnErr)
          g = - Aprod(r,2);
          if (~isempty(c))
             g = g + c;
          end
          if (mu > 0)
             g = g + mu * x;
          end

          s    = x - xOld;
          y    = g - gOld;
          sts  = s'*s;
          sty  = s'*y;
          if   sty <= 0,  bbScaling = stepMax;
          else            bbScaling = min( stepMax, max(stepMin, sts/sty) );
          end
       else
          bbScaling = min( stepMax, bbScaling );
       end


       %---------------------------------------------------------------
       % Update a Hessian approximation whenever possible
       %---------------------------------------------------------------
       if ((qnType ~= QN_NONE) && (~lnErr))

          % Step 1. Determine support
          flagSupport = 0;
          if ((abs(xNorm1 - tau) / max(1,tau)) > 1e-8)
             support = false(n,1);    % Interior of the L1 ball
          else
             support = (absx > 1e-9);
          end

          % Step 2. Check if the support has changed
          if (any(support ~= supportOld))
             % Support has changed
             supportChange = iter;

             if (all(support(supportOld)) && all(sign(x(supportOld)) == xOld(supportOld)))
                flagSupport = 2; % Changed to parent face
             else
                flagSupport = 1;
             end
          else
             % Same support, different signs?
             if (any(support) && any(sign(x(support)) ~= sign(xOld(support))))
                flagSupport = 1; % Changed
                supportChange = iter;
             end
          end

          % Step 3. Check if support change type allows Hessian update
          if (iter > 1)
             if (flagSupport ~= 0)
                H = []; % Reset Hessian whenever the support changes
             end

             switch qnType
                case QN_CONSECUTIVE
                   flagUpdateHessian = (flagSupport == 0);
                   flagUseHessian    = flagUpdateHessian;

                case QN_PARENT_FAST
                   flagUpdateHessian = ((flagSupport == 0) || (flagSupport == 2));
                   flagUseHessian    = flagUpdateHessian;

                case QN_PARENT_SLOW
                   flagUpdateHessian = ((flagSupport == 0) || (flagSupport == 2));
                   flagUseHessian    = (flagSupport == 0);

                otherwise
                   flagUpdateHessian = false;
                   flagUseHessian    = false;
             end

          else
             flagUpdateHessian = false;
             flagUseHessian    = false;
          end

          % Step 4. Check if face type allows Hessian update
          if (flagUpdateHessian)
             if (~any(support))
                if (qnAllowFull == 0)
                   % Disable Hessian update
                   flagUpdateHessian = false;
                end
             else
                nSupport = sum(support); % Also used in step 6.
                if ((nSupport == 1) || ((qnStrict) && (nSupport > m)))
                   % Disable Hessian update for vertices and high-dim. faces
                   flagUpdateHessian = false;
                end
             end
          end

          % Step 5. Check self-projection condition
          if (flagUpdateHessian)
             if (any(support)) % Boundary of crosspolytope
                d = -g;
                s1 = ((d <  0) & (x > 0)) | ((d >  0) & (x < 0));
                s2 = ((d <= 0) & (x < 0)) | ((d >= 0) & (x > 0));
                s3 = ~(s1 | s2); % Zero entries in x

                % Compute the absolute sum of the sets
                sum1 = sum(abs(d(s1))); % Decreases the one norm
                sum2 = sum(abs(d(s2))); % Increases the one norm
                sum3 = sum(abs(d(s3))); % Increases the one norm

                if (sum1 > sum2+sum3)
                   % Move into the polytope
                   flagSelfProj = false;
                elseif (sum1 < sum2+sum3)
                   % Move away from the boundary
                   flagSelfProj = (max(abs(d(s3))) <= (sum1+sum2) / sum(support));
                else
                   % Move across the face
                   flagSelfProj = (sum3 == 0);
                end
             else % Interior of crosspolytope
                flagSelfProj = true;
             end

             if (~flagSelfProj)
                flagUpdateHessian = false;
                statusMsg = [statusMsg,'-'];
             end
          end

          % Step 6. Initialize and/or update Hessian approximation
          if (flagUpdateHessian)
             if (isempty(H))
                if (~any(support))
                   % Interior of crosspolytope
                   H = lbfgsinit(n,lbfgsHist,1e-3);
                else
                   % Boundary of crosspolytope
                   signs = sign(x(support));
                   H = lbfgsinit(nSupport-1,lbfgsHist,1e-3);
                end
             end

             % Update Hessian approximation
             if (~any(support))
                H = lbfgsupdate(H, 1, s, gOld, g);
             else
                if (isempty(sqrt1))
                   % Compute all square root factors
                   sqrt1 = sqrt(1 ./ ((1:n) .* (2:n+1)))'; % sqrt(1 / (i*(i+1)))
                   sqrt2 = sqrt((1:n) ./ (2:n+1))';        % sqrt(i / (i + 1.0))
                end

                H = lbfgsupdate(H, 1, ...
                   productBMex(signs.*s(support)   ,1,sqrt1,sqrt2), ...
                   productBMex(signs.*gOld(support),1,sqrt1,sqrt2), ...
                   productBMex(signs.*g(support)   ,1,sqrt1,sqrt2));
             end
          else
             flagUseHessian = false;
             H = [];
          end
       end % qnType ~= QN_NONE


    catch err % Detect matrix-vector multiply limit error
       if strcmp(err.identifier,'SPGL1:MaximumMatvec')
         stat = EXIT_MATVEC_LIMIT;
         iter = iter - 1;
         x = xOld;  f = fOld;  g = gOld;  r = rOld;
         break;
       else
         rethrow(err);
       end
    end


    %------------------------------------------------------------------
    % Update function history.
    %------------------------------------------------------------------
    if singleTau || f > sigma^2 / 2 % Don't update if superoptimal.
       lastFv(mod(iter,nPrevVals)+1) = f;
       if fBest > f
          fBest = f;
          xBest = x;
       end
    end

end % while 1


% Restore best solution (only if solving single problem).
if singleTau && f > fBest
   rNorm = sqrt(2*fBest);
   printf('\n Restoring best iterate to objective %13.7e\n',rNorm);
   x = xBest;
   r = b - Aprod(x,1);
   g =   - Aprod(r,2);
   if (~isempty(c))
      g = g + c;
   end
   if (mu > 0)
      g = g + mu * x;
   end


   gNorm = norm(g,inf);
   rNorm = norm(r,  2);
end

% Final cleanup before exit.
info.tau           = tau;
info.f             = f;
info.fBest         = fBest;
info.fDual         = fDual;
info.rNorm         = rNorm;
info.rGap          = rGap;
info.gNorm         = gNorm;
info.rGap          = rGap;
info.stat          = stat;
info.iter          = iter;
info.iterH         = iterH;
info.nProdA        = nProdA;
info.nProdAt       = nProdAt;
info.nNewton       = nNewton;
info.timeProject   = timeProject;
info.timeMatProd   = timeMatProd;
info.options       = options;
info.timeTotal     = toc;
info.supportChange = supportChange; % Iteration of last support change

if (history)
   info.historyTime = historyTime(1:iter+1);
   info.historyFun  = historyFun(1:iter+1);
   info.historyDual = historyDual(1:iter+1);
   info.historyGap  = historyGap(1:iter+1);
   info.historyQN   = historyQN(1:iter+1);
end

% Print final output.
switch (stat)
   case EXIT_OPTIMAL
      printf('\n EXIT -- Optimal solution found\n')
   case EXIT_ITERATIONS
      printf('\n ERROR EXIT -- Too many iterations\n');
   case EXIT_ROOT_FOUND
      printf('\n EXIT -- Found a root\n');
   case {EXIT_BPSOL_FOUND}
      printf('\n EXIT -- Found a BP solution\n');
   case {EXIT_LEAST_SQUARES}
      printf('\n EXIT -- Found a least-squares solution\n');
   case EXIT_LINE_ERROR
      printf('\n ERROR EXIT -- Linesearch error (%i)\n',lnErr);
   case EXIT_SUBOPTIMAL_BP
      printf('\n EXIT -- Found a suboptimal BP solution\n');
   case EXIT_MATVEC_LIMIT
      printf('\n EXIT -- Maximum matrix-vector operations reached\n');
   otherwise
      error('Unknown termination condition\n');
end
printf('\n');
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Products with A',nProdA,'','Total time   (secs)',info.timeTotal);
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Products with A''',nProdAt,'','Project time (secs)',timeProject);
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Newton iterations',nNewton,'','Mat-vec time (secs)',timeMatProd);
printf(' %-20s:  %6i\n', ...
   'Line search its',nLineTot);
printf('\n');



% =================================================================== %
% NESTED FUNCTIONS.  These share some vars with workspace above.      %
% =================================================================== %

function z = Aprod(x,mode)
   if (nProdA + nProdAt >= maxMatvec)
     error('SPGL1:MaximumMatvec','');
   end

   tStart = toc;
   if mode == 1
      nProdA = nProdA + 1;
      if   explicit, z = A*x;
      else           z = A(x,1);
      end
   elseif mode == 2
      nProdAt = nProdAt + 1;
      if   explicit, z = A'*x;
      else           z = A(x,2);
      end
   else
      error('Wrong mode!');
   end
   timeMatProd = timeMatProd + (toc - tStart);
end % function Aprod

% =================================================================== %

function printf(varargin)
  if logLevel > 0
     fprintf(fid,varargin{:});
  end
end % function printf


function x = project(x, tau, weights)
   tStart = toc;
   x = oneProjector(x,weights,tau);
   timeProject = timeProject + (toc - tStart);
end % function project


function v = option(field,value)
   if (isfield(options,field))
      v = options.(field);
   else
      v = value;
   end
end


% =================================================================== %
% End of nested functions.                                            %
% =================================================================== %

end % function spg

% =================================================================== %
% PRIVATE FUNCTIONS.                                                  %
% =================================================================== %

function [fNew,xNew,rNew,iter,err] = spgLine(f,x,d,gtd,fMax,Aprod,b,c,mu)
% Nonmonotone linesearch.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;
maxIts = 10;
step   = 1;
iter   = 0;
gamma  = 1e-4;
gtd    = real(gtd); % Need to use real when dealing with complex numbers
                    % TODO: Do not use -abs(.)
while 1
    % Evaluate trial point and function value.
    xNew = x + step*d;
    rNew = b - Aprod(xNew,1);
    fNew = rNew'*rNew / 2;
    if (~isempty(c))
       fNew = fNew + c'*xNew;
    end
    if (mu > 0)
       fNew = fNew + (mu/2) * (xNew'*xNew);
    end


    % Check exit conditions.
    if fNew < fMax + gamma*step*gtd  % Sufficient descent condition.
       err = EXIT_CONVERGED;
       break
    elseif  iter >= maxIts           % Too many linesearch iterations.
       err = EXIT_ITERATIONS;
       break
    end

    % New linesearch iteration.
    iter = iter + 1;

    % Safeguarded quadratic interpolation.
    if step <= 0.1
       step  = step / 2;
    else
       tmp = (-gtd*step^2) / (2*(fNew-f-step*gtd));
       if tmp < 0.1 || tmp > 0.9*step || isnan(tmp)
          tmp = step / 2;
       end
       step = tmp;
    end

end % while 1

end % function spgLine


% =================================================================== %


function [fNew,xNew,rNew,iter,step,err] = ...
    spgLineCurvy(x,g,fMax,Aprod,b,c,mu,project,tau,weights)
% Projected backtracking linesearch.
% On entry,
% g  is the (possibly scaled) steepest descent direction.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;
EXIT_NODESCENT  = 2;
gamma  = 1e-4;
maxIts = 10;
step   =  1;
sNorm  =  0;
scale  =  1;      % Safeguard scaling.  (See below.)
nSafe  =  0;      % No. of safeguarding steps.
iter   =  0;
debug  =  false;  % Set to true to enable log.
n      =  length(x);

if debug
   fprintf(' %5s  %13s  %13s  %13s  %8s\n',...
           'LSits','fNew','step','gts','scale');
end

while 1

    % Evaluate trial point and function value.
    xNew     = project(x - step*scale*g, tau, weights);
    rNew     = b - Aprod(xNew,1);
    fNew     = rNew'*rNew / 2;
    if (~isempty(c))
       fNew = fNew + c'*xNew;
    end
    if (mu > 0)
       fNew = fNew + (mu/2) * (xNew'*xNew);
    end

    s   = xNew - x;
    gts = scale * real(g' * s);
    if gts >= 0
       err = EXIT_NODESCENT;
       break
    end

    if debug
       fprintf(' LS %2i  %13.7e  %13.7e  %13.6e  %8.1e\n',...
               iter,fNew,step,gts,scale);
    end

    % 03 Aug 07: If gts is complex, then should be looking at -abs(gts).
    % 13 Jul 11: It's enough to use real part of g's (see above).
    if fNew < fMax + gamma*step*gts
%   if fNew < fMax - gamma*step*abs(gts)  % Sufficient descent condition.
       err = EXIT_CONVERGED;
       break
    elseif iter >= maxIts                 % Too many linesearch iterations.
       err = EXIT_ITERATIONS;
       break
    end

    % New linesearch iteration.
    iter = iter + 1;
    step = step / 2;

    % Safeguard: If stepMax is huge, then even damped search
    % directions can give exactly the same point after projection.  If
    % we observe this in adjacent iterations, we drastically damp the
    % next search direction.
    % 31 May 07: Damp consecutive safeguarding steps.
    sNormOld  = sNorm;
    sNorm     = norm(s) / sqrt(n);
    %   if sNorm >= sNormOld
    if abs(sNorm - sNormOld) <= 1e-6 * sNorm
       gNorm = norm(g) / sqrt(n);
       scale = sNorm / gNorm / (2^nSafe);
       nSafe = nSafe + 1;
    end

end % while 1

end % function spgLineCurvy
