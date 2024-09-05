% Copyright contributors to the Hybrid-BPDN-solver project

function table_correlated()

% Run SPGL1 on difficult problems, plot info.xNorm1
% How many of the problems reach 20 line-search iterations

% Correlated (Gaussian): 10,50 : 1:10, 91
% Correlated (Gaussian): 100   : 1:4,  91
% Problems 5,6, k=50: 10 runs

% Results from 9 and 10 show that low accuracy can also lead to
% much larger run times because of the large number of steps
% required to converge to the correct tau

% SPGL1 uses a heuristic, often works for simpler problems but
% can lead to a degradation in accuracy.

prefix = mfilepath(mfilename('fullpath'));
prefix_cache = [prefix,'cache',filesep];
if (~exist(prefix_cache, 'dir'))
   [status,msg,msgID] = mkdir(prefix_cache);
end


% ---------------------------------------------------------------------
% Collect all results
% ---------------------------------------------------------------------

%muVals = [0.8,0.85,0.9,0.95,0.98,0.99,0.995];
muVals = [0.9,0.95,0.98,0.99,0.995];
nTrials = 10;

methods         = [1:8,90,91]; % 92, 32 for correlated
problems        = 1:7;
sparsityVals    = [10,50,100];
sparsityIndices = {[1,2,3], ...
                   [1,2,3], ...
                   [1,2,3], ...
                   [2], ...
                   [2], ...
                   [2], ...
                   [2]};


filename = [prefix_cache, 'correlated_problems_combined.mat'];
if (~exist(filename,'file'))
   data = cell(length(sparsityVals), length(problems),length(methods));
   for problemIdx = 1:length(problems)
      fprintf('Loading data for problem %d . . .\n', problems(problemIdx));
      problem = problems(problemIdx);

      for sparsityIdx = sparsityIndices{problemIdx}
         sparsity = sparsityVals(sparsityIdx);

         for methodIdx = 1:length(methods)
            method = methods(methodIdx);

            info = getCorrelatedData(sparsity, problem, method, nTrials, muVals);

            data{sparsityIdx, problemIdx, methodIdx} = info;
         end
      end
   end
   save(filename,'data');
else
   data = load(filename);
   data = data.data;
end


% ---------------------------------------------------------------------
% TABLE 1 -  Root-finding time and speedup
% ---------------------------------------------------------------------
methodVals = [5,6,1,2,7,8,3,4,91];

% Prepare the information per block
blockInfo = {};
blockLabel = {};

% Blocks 1-3: Noiseless, runtimes averaged over the three support types
for sparsityIdx = 1:3
   runtime = zeros(length(muVals), length(methodVals));
   for i=1:length(methodVals)
      methodIdx = find(methods == methodVals(i), 1, 'first');
      for problemIdx = 1:3
         info = data{sparsityIdx, problemIdx, methodIdx};
         runtime(:,i) = runtime(:,i) + mean(info.runtime,1)';
      end
   end
   runtime = runtime / 3; % Average over three methods
   blockInfo{end+1} = runtime;
   blockLabel{end+1} = sprintf('Sparsity $k=%d$',sparsityVals(sparsityIdx));
end

% Blocks 4-7: Noisy
for block=1:4
   noiseLevels = [1,2,5,10];
   problemIdx  = block + 3; % Problems 4 to 7
   sparsityIdx = 2;         % Sparsity level 50

   if (block == 2), continue; end % Skip 2%

   runtime = zeros(length(muVals), length(methodVals));
   for i=1:length(methodVals)
      methodIdx = find(methods == methodVals(i), 1, 'first');
      info = data{sparsityIdx, problemIdx, methodIdx};
      runtime(:,i) = mean(info.runtime,1)';
   end
   blockInfo{end+1} = runtime;
   blockLabel{end+1} = sprintf('Sparsity $k=50$, noise level %d\\%%', noiseLevels(block));
end

% Prepare table output stream
fp = 1;

fprintf(fp,'\\begin{tabular}{rcrcrrrcrrrcrrrcrrr}\n');
fprintf(fp,['$\\gamma$ && SPGL1 && \\multicolumn{3}{c}{Root finding}',...
                               '&& \\multicolumn{3}{c}{Root finding}',...
                               '&& \\multicolumn{3}{c}{Lasso}',...
                               '&& \\multicolumn{3}{c}{Lasso}\\\\\n']);
fprintf(fp,['&& runtime && \\multicolumn{3}{c}{tolerance $10^{-4}$}',...
                       '&& \\multicolumn{3}{c}{tolerance $10^{-6}$}',...
                       '&& \\multicolumn{3}{c}{tolerance $10^{-4}$}',...
                       '&& \\multicolumn{3}{c}{tolerance $10^{-6}$}\\\\\n']);
fprintf(fp,'\\cline{1-1}\\cline{3-3}\\cline{5-7}\\cline{9-11}\\cline{13-15}\\cline{17-19}\n');


for blockIdx=1:length(blockInfo)
   runtime = blockInfo{blockIdx};

   for muIdx = 1:length(muVals)
      str = sprintf('%.3f', 1-muVals(muIdx));
      fprintf(fp,'%5s',str(2:end));

      v = runtime(muIdx,:);

      % SPGL1
      fprintf(fp,'&& %5.1f', v(end));

      % Other methods
      for i=1:2:length(methodVals)-1
         fprintf(fp,'&& %5.1f & %5.1f & %4d', v(i),v(i+1),round(100*(v(i)-v(i+1)) / v(i)));
      end

      if (muIdx == length(muVals))
         str = '[4pt]';
      else
         str = '';
      end
      fprintf(fp,'\\\\%s\n', str);
   end

   if (blockIdx ~= length(blockInfo))
      str = '[8pt]';
   else
      str = '';
   end
   fprintf(fp,'\\multicolumn{19}{c}{({\\bf{%s}}) %s}\\\\%s\n',...
      char('a' + blockIdx-1),blockLabel{blockIdx},str);
end

% Table footer
fprintf(fp,'\\end{tabular}');

% Close table output stream if needed
fprintf(fp,'\n');


% ---------------------------------------------------------------------
% TABLE 2 -  Status
% ---------------------------------------------------------------------
methodVals = [5,6,1,2,7,8,3,4,91];

% Prepare the information per block
blockInfo = {};
blockLabel = {};

% Blocks 1-3: Noiseless, runtimes averaged over the three support types
for sparsityIdx = 1:3
   status = zeros(length(muVals), length(methodVals));
   for i=1:length(methodVals)
      methodIdx = find(methods == methodVals(i), 1, 'first');
      for problemIdx = 1:3
         info = data{sparsityIdx, problemIdx, methodIdx};
         status(:,i) = status(:,i) + sum(info.status(:,:) ~= 4)';
      end
   end
   %status = status / 30; % Average over three methods
   blockInfo{end+1} = status;
   blockLabel{end+1} = sprintf('Sparsity $k=%d$',sparsityVals(sparsityIdx));
end

% Blocks 4-7: Noisy
for block=1:4
   noiseLevels = [1,2,5,10];
   problemIdx  = block + 3; % Problems 4 to 7
   sparsityIdx = 2;         % Sparsity level 50

   if (block == 2), continue; end % Skip 2%

   status = zeros(length(muVals), length(methodVals));
   for i=1:length(methodVals)
      methodIdx = find(methods == methodVals(i), 1, 'first');
      info = data{sparsityIdx, problemIdx, methodIdx};
      status(:,i) = status(:,i) + sum(info.status(:,end) ~= 4)';
   end
   blockInfo{end+1} = status;
   blockLabel{end+1} = sprintf('Sparsity $k=50$, noise level %d\\%%', noiseLevels(block));
end

% Prepare table output stream
fp = 1;

fprintf(fp,'\\begin{tabular}{rcrcrrcrrcrrcrr}\n');
fprintf(fp,['$\\gamma$ && SPGL1 && \\multicolumn{2}{c}{Root finding}',...
                               '&& \\multicolumn{2}{c}{Root finding}',...
                               '&& \\multicolumn{2}{c}{Lasso}',...
                               '&& \\multicolumn{2}{c}{Lasso}\\\\\n']);
fprintf(fp,['&& error && \\multicolumn{2}{c}{tol. $10^{-4}$}',...
                       '&& \\multicolumn{2}{c}{tol. $10^{-6}$}',...
                       '&& \\multicolumn{2}{c}{tol. $10^{-4}$}',...
                       '&& \\multicolumn{2}{c}{tol. $10^{-6}$}\\\\\n']);
fprintf(fp,'\\cline{1-1}\\cline{3-3}\\cline{5-6}\\cline{8-9}\\cline{11-12}\\cline{14-15}\n');


for blockIdx=1:length(blockInfo)
   runtime = blockInfo{blockIdx};

   for muIdx = 1:length(muVals)
      str = sprintf('%.3f', 1-muVals(muIdx));
      fprintf(fp,'%5s',str(2:end));

      v = runtime(muIdx,:);

      % SPGL1
      fprintf(fp,'&& %d', v(end));

      % Other methods
      for i=1:2:length(methodVals)-1
         fprintf(fp,'&& %d & %d', v(i),v(i+1));
      end

      if (muIdx == length(muVals))
         str = '[4pt]';
      else
         str = '';
      end
      fprintf(fp,'\\\\%s\n', str);
   end

   if (blockIdx ~= length(blockInfo))
      str = '[8pt]';
   else
      str = '';
   end
   fprintf(fp,'\\multicolumn{15}{c}{({\\bf{%s}}) %s}\\\\%s\n',...
      char('a' + blockIdx-1),blockLabel{blockIdx},str);
end

% Table footer
fprintf(fp,'\\end{tabular}');

% Close table output stream if needed
fprintf(fp,'\n');



% ---------------------------------------------------------------------
% TABLE 3 -  Median optimality
% ---------------------------------------------------------------------

methodVals = [5,6,1,2,7,8,3,4,91];

% Prepare the information per block
blockInfo = {};
blockLabel = {};

% Blocks 1-3: Noiseless, runtimes averaged over the three support types
for sparsityIdx = 1:3
   optimality = zeros(length(muVals), length(methodVals));
   for i=1:length(methodVals)
      methodIdx = find(methods == methodVals(i), 1, 'first');

      % Problem indices 1 through 3
      info1 = data{sparsityIdx, 1, methodIdx};
      info2 = data{sparsityIdx, 2, methodIdx};
      info3 = data{sparsityIdx, 3, methodIdx};

      %optimality(:,i) = median([info1.relgap; info2.relgap; info3.relgap],1)';
      optimality(:,i) = exp(mean(log([info1.relgap; info2.relgap; info3.relgap]),1))';
      %optimality(:,i) = max([info1.relgap; info2.relgap; info3.relgap],[],1)';

   end
   %status = status / 30; % Average over three methods
   blockInfo{end+1} = optimality;
   blockLabel{end+1} = sprintf('Sparsity $k=%d$',sparsityVals(sparsityIdx));
end

% Blocks 4-7: Noisy
for block=1:4
   noiseLevels = [1,2,5,10];
   problemIdx  = block + 3; % Problems 4 to 7
   sparsityIdx = 2;         % Sparsity level 50

   if (block == 2), continue; end % Skip 2%

   optimality = zeros(length(muVals), length(methodVals));
   for i=1:length(methodVals)
      methodIdx = find(methods == methodVals(i), 1, 'first');
      info = data{sparsityIdx, problemIdx, methodIdx};
      %optimality(:,i) = median(info.relgap,1)';
      optimality(:,i) = exp(mean(log(info.relgap),1))';
      %optimality(:,i) = max(info.relgap,[],1)';
   end
   blockInfo{end+1} = optimality;
   blockLabel{end+1} = sprintf('Sparsity $k=50$, noise level %d\\%%', noiseLevels(block));
end


% Prepare table output stream
fp = 1;

fprintf(fp,'\\begin{tabular}{rcrcrrcrrcrrcrr}\n');
fprintf(fp,['$\\gamma$ && SPGL1 && \\multicolumn{2}{c}{Root finding}',...
                               '&& \\multicolumn{2}{c}{Root finding}',...
                               '&& \\multicolumn{2}{c}{Lasso}',...
                               '&& \\multicolumn{2}{c}{Lasso}\\\\\n']);
fprintf(fp,['&& rel. gap && \\multicolumn{2}{c}{tol. $10^{-4}$}',...
                       '&& \\multicolumn{2}{c}{tol. $10^{-6}$}',...
                       '&& \\multicolumn{2}{c}{tol. $10^{-4}$}',...
                       '&& \\multicolumn{2}{c}{tol. $10^{-6}$}\\\\\n']);
fprintf(fp,'\\cline{1-1}\\cline{3-3}\\cline{5-6}\\cline{8-9}\\cline{11-12}\\cline{14-15}\n');


for blockIdx=1:length(blockInfo)
   relgap = blockInfo{blockIdx};

   for muIdx = 1:length(muVals)
      str = sprintf('%.3f', 1-muVals(muIdx));
      fprintf(fp,'%5s',str(2:end));

      v = relgap(muIdx,:);

      % SPGL1
      str = sprintf('%5.1e', v(end)); str = [str(1:5),str(7:end)];
      fprintf(fp,'&& %s', str);

      % Other methods
      for i=1:2:length(methodVals)-1
         str1 = sprintf('%5.1e', v(i  )); str1 = [str1(1:5),str1(7:end)];
         str2 = sprintf('%5.1e', v(i+1)); str2 = [str2(1:5),str2(7:end)];
         fprintf(fp,'&& %s & %s', str1, str2);
      end

      if (muIdx == length(muVals))
         str = '[4pt]';
      else
         str = '';
      end
      fprintf(fp,'\\\\%s\n', str);
   end

   if (blockIdx ~= length(blockInfo))
      str = '[8pt]';
   else
      str = '';
   end
   fprintf(fp,'\\multicolumn{15}{c}{({\\bf{%s}}) %s}\\\\%s\n',...
      char('a' + blockIdx-1),blockLabel{blockIdx},str);
end

% Table footer
fprintf(fp,'\\end{tabular}');

% Close table output stream if needed
fprintf(fp,'\n');


% ---------------------------------------------------------------------
% FIGURE 1 -  Runtime
% ---------------------------------------------------------------------

%methodVals = [5,6,1,2,7,8,3,4,91];
methodVals  = [91,1,2,5,6,3,4,7,8];
compareVals = [ 0,0,1,0,5,0,3,0,7];
typeVals    = [1,1,1,1,1,2,2,2,2];
typeName    = {'BP$_{\sigma}$','LS$_{\tau}$'};
acc         = [1e-6, 1e-6,1e-6,1e-4,1e-4,1e-6,1e-6,1e-4,1e-4];
methodName  = {'\spgl', ...
               'Original','Hybrid','Original','Hybrid', ...
               'Original','Hybrid','Original','Hybrid'};
bounds      = [-inf, 1e-6, 1e-5, 1e-4,1e-3,1e-2,inf];

% Problem indices refer to the on-support distribution, sparsity values
% refers to the sparsity of x0. For problems 1-3 we have three sparsity
% levels, for problems 4-7 we have one sparsity level each.
nElements  =             3 * 3 * length(muVals) * nTrials;
nElements  = nElements + 3 * 1 * length(muVals) * nTrials;
runtime    = zeros(nElements, length(methodVals));
status     = zeros(nElements, length(methodVals));
optimality = zeros(nElements, length(methodVals));

for i=1:length(methodVals)
   methodIdx = find(methods == methodVals(i), 1, 'first');

   offset = 0;
   for problemIdx = 1:7
      if (problemIdx < 4)
         sparsityValues = 1:3;
      else
         sparsityValues = 2;
      end

      for sparsityIdx = sparsityValues;
         info = data{sparsityIdx, problemIdx, methodIdx};
         runtime(offset+(1:nTrials*length(muVals)),i) = info.runtime(:);
         status(offset+(1:nTrials*length(muVals)),i) = info.status(:);
         optimality(offset+(1:nTrials*length(muVals)),i) = info.relgap(:);
         offset = offset + nTrials * length(muVals);
      end
   end
end


% Prepare table output stream
fp = 1;

% Table header
fprintf(fp,'\\begin{tabular}{lllr%s}\n',repmat('r',1,length(bounds)-1));
fprintf(fp,'\\hline\n');
fprintf(fp,'&&&&\\multicolumn{%d}{c}{Relative duality gap}\\\\\n',length(bounds)-1);
fprintf(fp,'\\cline{5-%d}\n',4+length(bounds)-1);
fprintf(fp,'\\\\[-10pt]\n');
fprintf(fp,'Type & Method & Tol. & \\multicolumn{1}{l}{Time} ');
for i=1:length(bounds)-1
   if (isinf(bounds(i)))
      fprintf(fp,'& $\\leq 10^{%d}$',log10(bounds(i+1)));
   elseif (isinf(bounds(i+1)))
      fprintf(fp,'& $> 10^{%d}$',log10(bounds(i)));
   else
      fprintf(fp,'& $10^{(%d,%d]}$',log10(bounds(i)),log10(bounds(i+1)));
   end
end
fprintf(fp,'\\\\\n');
fprintf(fp,'\\hline\n');


type = 0;
for i=1:length(methodVals)
   % Type
   if (type == typeVals(i))
      fprintf(fp,'%-14s','');
   else
      type = typeVals(i);
      fprintf(fp,'%-14s', typeName{type});
   end

   % Method
   fprintf(fp,'& %-10s', methodName{i});

   % Optimality tolerance
   fprintf(fp,'& $10^{%d}$', log10(acc(i)));

   % Runtime
   t = sum(runtime(:,i));

   h = floor(t/3600);
   m = floor((t - 3600*h)/60);
   s = round(t - 3600*h - 60*m);

   fprintf(fp,'& %2dh%02d', h,m);


   for j=1:length(bounds)-1
      s = 100*sum((optimality(:,i) > bounds(j)) & (optimality(:,i) <= bounds(j+1)))/size(optimality,1);

      if (s == 0)
         fprintf(fp,'& --');
      elseif (s < 10)
         fprintf(fp,'& %5.1f', s);
      else
         fprintf(fp,'& %5.0f', s);
      end
   end

   fprintf(fp,'\\\\\n');
end

% Table footer
fprintf(fp,'\\hline\n');
fprintf(fp,'\\end{tabular}');

% Close table output stream if needed
fprintf(fp,'\n');



for i=1:length(methodVals)
   idx = find(methodVals==compareVals(i));
   if (~isempty(idx))
      t1 = sum(runtime(:,i));
      t2 = sum(runtime(:,idx));
      fprintf('Runtime improvement %2d%%\n', round(100*((t2-t1) / t2)));
   end
end

return



labels = {'BP$_{\sigma}$ ($10^{-4}$)', ...
          'BP$_{\sigma}$ ($10^{-6}$)', ...
          'LS$_{\tau}$ ($10^{-4}$)', ...
          'LS$_{\tau}$ ($10^{-6}$)'};

fp = 1;
fprintf(fp,'\\begin{tabular}{lrrr}\n');
fprintf(fp,'Setting & Total time & Total time & Reduction (\\%%)\\\\\n');
fprintf(fp,'original & hybrid, & \\\\\n');
for j=1:4
   i = 2*j-1;
   t1 = sum(runtime(:,i));
   t2 = sum(runtime(:,i+1));

   h1 = floor(t1/3600);
   m1 = floor((t1 - 3600*h1)/60);
   s1 = round(t1 - 3600*h1 - 60*m1);

   h2 = floor(t2/3600);
   m2 = floor((t2 - 3600*h2)/60);
   s2 = round(t2 - 3600*h2 - 60*m2);

   fprintf(fp,'%s & %2d:%02d:%02d & %2d:%02d:%02d &%.0f\\%%\\\\\n',...
           labels{j}, h1,m1,s1, h2,m2,s2, 100*(t1-t2)./t1);
end
fprintf(fp,'\\end{tabular}');


fprintf(fp,'\n');
for j=1:9
   fprintf('%d   ',j);
   for i = 1:length(bound)-1
      %fprintf('  %5.1f', 100*sum((optimality(:,j) > bound(i)) & (optimality(:,j) <= bound(i+1)))/size(optimality,1));
      fprintf('  %5.1f', 100*sum((optimality(:,j) > bound(i)))/size(optimality,1));
   end
   fprintf('\n');
end



function info = getCorrelatedData(sparsity, problem, method, nTrials, muVals)

runtime    = zeros(nTrials,length(muVals));
xratio     = zeros(nTrials,length(muVals));
iterations = zeros(nTrials,length(muVals));
matprod    = zeros(nTrials,length(muVals));
nNewton    = zeros(nTrials,length(muVals));
tau        = zeros(nTrials,length(muVals));
status     = zeros(nTrials,length(muVals));
rGap       = zeros(nTrials,length(muVals));

for i=1:nTrials
   for muIdx=1:length(muVals)
      mu = muVals(muIdx);

      %fprintf('k = %2d, problem = %d/3, method = %2d, mu = %f, idx = %d\n', sparsity,problem,method,mu,i);
      d = run_problem_correlated(problem,method,sparsity,mu,i);

      if (iscell(d.data.info))
         iter=zeros(1,length(d.data.info));
         for j=1:length(d.data.info)
            iter(j) = d.data.info{j}.iter;
         end
      end

      runtime(i,muIdx)    = d.runtime;
      xratio(i,muIdx)     = d.xOffSupport / (d.xOffSupport + d.xOnSupport);
      iterations(i,muIdx) = d.iterations;
      matprod(i,muIdx)    = d.matprod;
      nNewton(i,muIdx)    = d.nNewton;
      tau(i,muIdx)        = d.tau;
      status(i,muIdx)     = d.status(end);
      rGap(i,muIdx)       = d.rGap(end);
   end
end

info = struct();
info.runtime    = runtime;
info.xratio     = xratio;
info.iterations = iterations;
info.matprod    = matprod;
info.nNewton    = nNewton;
info.tau        = tau;
info.status     = status;
info.relgap     = rGap;
