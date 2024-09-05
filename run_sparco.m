% Copyright contributors to the Hybrid-BPDN-solver project

function data = run_sparco(problemIdx, scale, method, optTol)

% Determine the filename for caching
prefix = mfilepath(mfilename('fullpath'));
prefix_cache = [prefix,'cache',filesep,'sparco',filesep];
if (~exist(prefix_cache, 'dir'))
   [status,msg,msgID] = mkdir(prefix_cache);
end

multiplier = 1;
if ((problemIdx == 603) || (problemIdx == 903) || (problemIdx == 701) || (problemIdx == 702))
   multiplier = 100;
end

scaleStr = sprintf('%.4f',scale);
scaleStr(2) = '_';
optTolStr = sprintf('1e%d',round(log10(optTol)));
filename = sprintf('Sparco-%03d-%s-m%d-%s.mat', problemIdx, scaleStr, method, optTolStr);
filename = [prefix_cache, filename];

if (exist(filename,'file'))
   d = load(filename);
   data = d.data;
else
   try
      p = generateProblem(problemIdx);
      p.b = p.b * multiplier;
   catch exception
      warning('Code requires external dependency Sparco.')
      rethrow(exception)
   end
   sigma = scale * norm(p.b,2);


   % Method values
   % 1 - spgl1
   % 2 - spgl1, optTol = 1e-6, decTol = 1e-6
   % 3 - QN=0
   % 4 - QN=1

   if (method == 1)
      % SPGL1
      options = struct();
      options.optTol    = 1e-4;
      options.decTol    = 1e-4;
      options.verbosity = 0;
      [x,r,g,info] = spgl1(p.A,p.b,[],sigma,[],options);
      runtime = info.timeTotal;
      rGap    = info.rGap;

   elseif (method == 2)
      % SPGL1 - Tight
      options = struct();
      options.optTol    = 1e-6;
      options.decTol    = 1e-6;
      options.verbosity = 0;
      [x,r,g,info] = spgl1(p.A,p.b,[],sigma,[],options);
      runtime = info.timeTotal;
      rGap    = info.rGap;

   elseif (method == 6)
      % SPGL1 - Super tight
      options = struct();
      options.optTol    = 1e-9;
      options.decTol    = 1e-9;
      options.verbosity = 0;
      [x,r,g,info] = spgl1(p.A,p.b,[],sigma,[],options);
      runtime = info.timeTotal;
      rGap    = info.rGap;

   elseif (method == 3)
      % QN = 0
      options = struct();
      options.qnType     = 0;
      options.optTol     = optTol;
      options.optTolMinF = 1;
      options.verbosity  = 0;

      [x,r,g,info] = solve_rootfinding_v05(p.A,p.b,sigma,options);
      runtime = sum(info.runtime);
      rGap    = info.info{end}.rGap;

   elseif (method == 4)
      % QN = 1
      options = struct();
      options.qnType     = 1;
      options.optTol     = optTol;
      options.optTolMinF = 1;
      options.verbosity  = 0;

      [x,r,g,info] = solve_rootfinding_v05(p.A,p.b,sigma,options);
      runtime = sum(info.runtime);
      rGap    = info.info{end}.rGap;

   else
      error('Invalid method index');
   end

   data = struct();
   data.info = info;
   data.runtime = runtime;
   data.xNorm1  = norm(x,1);
   data.rGap    = rGap;
   data.rNorm   = norm(r,2);
   data.obj     = (r'*r) / 2;
   save(filename,'data');
end
