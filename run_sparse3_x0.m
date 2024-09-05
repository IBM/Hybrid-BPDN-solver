% Copyright contributors to the Hybrid-BPDN-solver project

function data = run_sparse3_x0(m,n,k, type, method, nInstances)

prefix = mfilepath(mfilename('fullpath'));
prefix_cache = [prefix, 'cache', filesep, 'sparse', filesep];
if (~exist(prefix_cache, 'dir'))
   [status,msg,msgID] = mkdir(prefix_cache);
end

filename = [prefix_cache, sprintf('sparse3_x0_%d_%d_%d_%d_%d.mat', m,n,k,type, method)];
if (exist(filename,'file'))
   data = load(filename);
   startInstance = data.instance + 1;
   data = data.data;

   if (nInstances > size(data,1))
      data(nInstances,end) = 0;
   end
else
   % Data fields:
   % 1. iterations
   % 2. runtime
   % 3. exit code
   % 4. f
   % 5. fDual
   % 6. nProdA
   % 7. nProdAt
   % 8. norm x0
   % 9. norm ||x0 - x||_1
   % 10. norm ||x0 - x||_2
   % 11. norm b
   data = zeros(nInstances, 11);
   startInstance = 1;
end

optTol = 1e-6;

switch (method)
   case 1
      options = struct();
      options.qnType    = 0;
      options.lsType    = 1;
      options.optTol    = optTol;
      options.verbosity = 0;
      options.optTolMinF = 1e-3;
      options.iterations = 50000;
      %options.optTolRel = true; % make sure optimality tolerance is always relative

   case 2
      options = struct();
      options.qnType    = 1;
      options.lsType    = 1;
      options.optTol    = optTol;
      options.verbosity = 0;
      options.optTolMinF = 1e-3;
      options.iterations = 50000;
      %options.optTolRel = true; % make sure optimality tolerance is always relative

   case 3
      options = struct();
      options.qnType    = 0;
      options.lsType    = 7;
      options.optTol    = optTol;
      options.verbosity = 0;
      options.optTolMinF = 1e-3;
      options.iterations = 50000;
      %options.optTolRel = true; % make sure optimality tolerance is always relative

   otherwise
      error('Invalid choice of method.');
end


t0 = tic();
for i = startInstance : nInstances
   instance = i; % Avoid saving the loop variable itself
   fprintf('Working on (%d,%dx%d,%d) instance %d . . .\n', type, m,n,k, instance);

   [A,b,x0] = generate_problem_sparse_x0(m,n,k,type,instance);


   [x,r,g,info] = solver_v05( A, b, [], 0, 0.99*norm(x0,1), [], [], options);

   data(instance,:) = [info.iter, info.timeTotal, info.stat, ...
                       info.f, info.fDual, info.nProdA, info.nProdAt, ...
                       norm(x0,1),norm(x0-x,1),norm(x0-x,2), norm(b,2)];

   if ((instance == nInstances) || (toc(t0) > 60))
      save(filename,'data','instance');
      t0 = tic();
   end
end
