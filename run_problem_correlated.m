% Copyright contributors to the Hybrid-BPDN-solver project

function info = run_problem_correlated(problem,method,k,mu,index)
% Solves PROBLEM type (0 = Gaussian, 1 = Uniform, 2 = Sign) with the
% given METHOD (see below). The problem instances is given by INDEX,
% with column correlation MU and sparsity K.

% Problem dimensions
m = 1024; n = 2048;


switch(problem)
   case {1,2,3,4,5,6,7}
      cache_directory = sprintf('correlated_problem%02d', problem);
   otherwise
      error('Invalid problem index.');
end

% Determine the filename for caching
prefix = mfilepath(mfilename('fullpath'));
prefix_cache = [prefix,'cache',filesep,cache_directory,filesep];
if (~exist(prefix_cache, 'dir'))
   [status,msg,msgID] = mkdir(prefix_cache);
end

filename = [prefix_cache, sprintf('soln_%d_%d_%05d_%d.mat',k,method,round(10000*mu),index)];
if (exist(filename,'file'))
   % Load the data
   data = load(filename);
   info = data.info;
else
   % Generate the problem
   [A,b,x0,sigma] = generate_problem_correlated(problem,m,n,k,mu,index);

   % Set default solver
   solver = 'v05';
   optTolMinF = 1e-3;

   switch method
      case 1
         options = struct();
         options.qnType     = 0;
         options.optTol     = 1e-6;
         options.verbosity  = 0;
         options.optTolMinF = optTolMinF;
         options.lsType     = 1;
         options.iterations = 100000;

         tau = [];
      case 2
         options = struct();
         options.qnType     = 1;
         options.optTol     = 1e-6;
         options.verbosity  = 0;
         options.optTolMinF = optTolMinF;
         options.lsType     = 1;
         options.iterations = 100000;

         tau = [];
      case 3
         % Same as 1, but use tau
         options = struct();
         options.qnType     = 0;
         options.optTol     = 1e-6;
         options.verbosity  = 0;
         options.optTolMinF = optTolMinF;
         options.lsType     = 1;
         options.iterations = 100000;

         d = run_problem_correlated(problem,1,k,mu,index); % Result of method #1
         tau = d.data.tau;

      case 4
         % Same as 2, but use tau
         options = struct();
         options.qnType     = 1;
         options.optTol     = 1e-6;
         options.verbosity  = 0;
         options.optTolMinF = optTolMinF;
         options.lsType     = 1;
         options.iterations = 100000;

         d = run_problem_correlated(problem,1,k,mu,index); % Result of method #1
         tau = d.data.tau;

      case 5
         % Same as 1, optTol 1e-4
         options = struct();
         options.qnType     = 0;
         options.optTol     = 1e-4;
         options.verbosity  = 0;
         options.optTolMinF = optTolMinF;
         options.lsType     = 1;
         options.iterations = 100000;

         tau = [];
      case 6
         % Same as 2, optTol 1e-4
         options = struct();
         options.qnType     = 1;
         options.optTol     = 1e-4;
         options.verbosity  = 0;
         options.optTolMinF = optTolMinF;
         options.lsType     = 1;
         options.iterations = 100000;

         tau = [];
      case 7
         % Same as 3, optTol 1e-4
         options = struct();
         options.qnType     = 0;
         options.optTol     = 1e-4;
         options.verbosity  = 0;
         options.optTolMinF = optTolMinF;
         options.lsType     = 1;
         options.iterations = 100000;

         d = run_problem_correlated(problem,1,k,mu,index); % Result of method #1
         tau = d.data.tau;

      case 8
         % Same as 4, optTol 1e-4
         options = struct();
         options.qnType     = 1;
         options.optTol     = 1e-4;
         options.verbosity  = 0;
         options.optTolMinF = optTolMinF;
         options.lsType     = 1;
         options.iterations = 100000;

         d = run_problem_correlated(problem,1,k,mu,index); % Result of method #1
         tau = d.data.tau;

      case 90
         % SPGL1
         options = struct();
         options.iterations = 100000;
         options.verbosity  = 0;

         solver = 'spgl1';
         tau = [];

      case 91
         % SPGL1
         options = struct();
         options.iterations = 100000;
         options.verbosity  = 0;
         options.optTol = 1e-6;
         options.decTol = 1e-6;

         solver = 'spgl1';
         tau = [];
   end

   if (isempty(tau))
      if (isempty(sigma))
         sigma = 0.01 * norm(b,2);
      end

      if (strcmp(solver,'v05'))
         [x,r,g,data] = solve_rootfinding_v05(A,b,sigma,options);
      elseif (strcmp(solver,'v04'))
         [x,r,g,data] = solve_rootfinding_v04(A,b,sigma,options);
      elseif (strcmp(solver,'spgl1'))
         [x,r,g,info] = spgl1(A,b,[],sigma,[],options);
         info = rmfield(info,{'xNorm1','rNorm2','lambda'});
         data = struct();
         data.info = info;
      else
         [x,r,g,data] = solve_rootfinding_v04(A,b,0.01*norm(b,2),options);
      end
      data.mode = 'sigma';

   else
      if (solver == 'v05')
         [x,r,g,info] = solver_v05(A,b,[],0,tau,[],[],options);
      else
         [x,r,g,info] = solver_v04(A,b,[],0,tau,[],[],options);
      end
      data = struct();
      data.info = info;
      data.options = options;
      data.tau = tau;
      data.mode = 'tau';
   end

   info = struct();
   info.data        = data;
   info.solver      = solver;
   info.xOffSupport = sum(abs(x(x0==0)));
   info.xOnSupport  = sum(abs(x(x0~=0)));

   save(filename,'info');
end


% Add commonly used fields
if (strcmp(info.solver,'spgl1'))
   info.runtime    = info.data.info.timeTotal;
   info.iterations = info.data.info.iter;
   info.matprod    = info.data.info.nProdA + info.data.info.nProdAt;
   info.nNewton    = info.data.info.nNewton;
   info.tau        = info.data.info.tau;
   info.status     = info.data.info.stat;
   info.rGap       = info.data.info.rGap;
elseif (strcmp(info.data.mode,'sigma'))
   info.runtime    = sum(info.data.runtime);
   info.iterations = sum(info.data.iterations);
   info.nNewton    = length(info.data.iterations);
   info.matprod    = 0;
   info.status     = zeros(1,length(info.data.info)-1);
   info.rGap       = zeros(1,length(info.data.info)-1);

   for i=2:length(info.data.info)
      info.matprod = info.matprod + info.data.info{i}.nProdA;
      info.matprod = info.matprod + info.data.info{i}.nProdAt;
      info.status(i-1) = info.data.info{i}.stat;
      info.rGap(i-1)   = info.data.info{i}.rGap;
   end
   info.tau  = info.data.tau;
else
   info.runtime = info.data.info.timeTotal;
   info.iterations = info.data.info.iter;
   info.matprod = info.data.info.nProdA + info.data.info.nProdAt;
   info.nNewton = 0;
   info.tau     = info.data.tau;
   info.status  = info.data.info.stat;
   info.rGap    = info.data.info.rGap;
end
