% Copyright contributors to the Hybrid-BPDN-solver project

function result = run_problem_duality(idxMu,idxSparsity, idxTau, idxSigma, instances)
%function runtime = experiment_duality_2()
% Measure the effect on runtime of the old and new duality gap.
% - We always use the best dual value
% - No quasi-Newton
% - Line search type number 2

prefix = mfilepath(mfilename('fullpath'));
prefix_cache = [prefix,'cache',filesep,'duality',filesep];
if (~exist(prefix_cache, 'dir'))
   [status,msg,msgID] = mkdir(prefix_cache);
end

filename = sprintf('%d_%d_%d_%d.mat', idxMu, idxSparsity, idxTau, idxSigma);
filename = [prefix_cache, filename];

if (exist(filename,'file'))
   data = load(filename);
   result = data.result;

   startIndex = result.instances + 1;

   if (size(result.runtime,1) < instances)
      result.runtime(instances,:)    = 0;
      result.iterations(instances,:) = 0;
      result.stat(instances,:)       = 0;
      result.sparsity(instances,:)   = 0;
      result.iterH(instances,:)      = 0;
   end
else
   result = struct();
   result.runtime    = zeros(instances,4);
   result.iterations = zeros(instances,4);
   result.stat       = zeros(instances,4);
   result.f          = zeros(instances,4);
   result.fDual      = zeros(instances,4);
   result.sparsity   = zeros(instances,4,9);
   result.iterH      = zeros(instances,4);

   startIndex = 1;
end


t0 = tic();

muValues = [0.1, 0.01, 0.001, 0.0001];
mu = muValues(idxMu);

for j = startIndex:instances
   %disp(j);
   [A,b,tau] = generate_problem_duality(idxSparsity, idxTau, idxSigma, j);

   fprintf('Running idxMu=%d, idxSparsity=%d, idxTau=%d, idxSigma=%d, instance %d/%d . . .\n', ...
           idxMu, idxSparsity, idxTau, idxSigma, j,instances);

   % Run the solver
   dualType = [0,1,0,1];
   qnType   = [false,false,true,true];
   for i=1:4
      options = struct();
      options.qnType  = qnType(i);
      options.lsType  = 2;
      options.dualType = dualType(i);
      options.verbosity = 0;
      [x,r,g,info] = solver_v05( A, b, [], mu, tau, [], [], options );

      result.runtime(j,i)    = info.timeTotal;
      result.iterations(j,i) = info.iter;
      result.stat(j,i)       = info.stat;
      result.f(j,i)          = info.f;
      result.fDual(j,i)      = info.fDual;
      result.sparsity(j,i,:) = sparsity(x);
      result.iterH(j,i)       = info.iterH;
      result.instances       = j;
   end

   %{
   disp(result.stat(j,:));
   disp(result.iterations(j,:));
   disp(result.iterH(j,:));
   disp(result.runtime(j,:));
   fprintf('%.6e %.6e\n', min(result.f(j,:)),max(result.fDual(j,:)));
   %}

   if ((j == instances) || (toc(t0) >= 120))
      save(filename,'result');
      t0 = tic();
   end
end
