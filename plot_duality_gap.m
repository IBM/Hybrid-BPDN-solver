% Copyright contributors to the Hybrid-BPDN-solver project

function plot_duality_gap()

% Determine the filename for caching
prefix = mfilepath(mfilename('fullpath'));
prefix_cache = [prefix,'cache',filesep,'duality_gap1',filesep];
if (~exist(prefix_cache, 'dir'))
   [status,msg,msgID] = mkdir(prefix_cache);
end

settings = [1e-1, 1, 5, 2,  4; ...  # muIdx=1
            1e-2, 2, 5, 1,  1; ...  # muIdx=2 (1,6,2,3)
            1e-3, 7, 6, 2, 10; ...  # muIdx=3
            1e-4, 5, 6, 3,  7];

optTol     = 1e-4;
optTolMinF = 1e-3; %1e-3; % optTolMinF = 1


optTolValues     = [1e-4, 1e-4, 1e-6];
optTolMinFValues = [1,    1e-3, 1e-3];

for idx=1:size(settings,1)
   for i=1:length(optTolValues)
      optTol     = optTolValues(i);
      optTolMinF = optTolMinFValues(i);

      figure(1); clf;
      title(sprintf('Setup = %d, opt = %d', idx, i));
      drawnow;

      filename = [prefix_cache, sprintf('Setting%d_opt%d.mat',idx,i)];
      if (exist(filename,'file'))
         data = load(filename);
         infoBaseline = data.infoBaseline;
         info0        = data.info0;
         info1        = data.info1;
      else
         mu          = settings(idx,1);
         idxSparsity = settings(idx,2);
         idxTau      = settings(idx,3);
         idxSigma    = settings(idx,4);
         j           = settings(idx,5);

         [A,b,tau] = generate_problem_duality(idxSparsity, idxTau, idxSigma, j);

         options = struct();
         options.optTol     = 1e-9;
         options.optTolMinF = optTolMinF;
         options.verbosity  = 1;
         options.lsType     = 2; % Regular
         options.history    = false;
         options.qnType     = 0;
         options.dualType   = 1;
         options.iterations = 500000;
         [x,r,g,infoBaseline] = solver_v04( A, b, [], mu, tau, [], [], options);

         options = struct();
         options.optTol     = optTol;
         options.optTolMinF = optTolMinF;
         options.verbosity  = 0;
         options.lsType     = 2; % Regular
         options.qnType     = 0;
         options.iterations = 100000;
         options.history    = true;

         options.dualType = 0;
         [x,r,g,info0] = solver_v04( A, b, [], mu, tau, [], [], options);

         options.dualType = 1;
         [x,r,g,info1] = solver_v04( A, b, [], mu, tau, [], [], options);

         save(filename, 'infoBaseline','info0','info1');
      end


      fp0 = -cummax(-info0.historyFun);
      fd0 =  cummax( info0.historyDual);

      fp1 = -cummax(-info1.historyFun);
      fd1 =  cummax( info1.historyDual);

      fdBase = max([infoBaseline.fDual,max(fd0),max(fd1)]);

      %{
      figure(1);
      plot(1:length(fp0), (fp0 - fd0) ./ max(fp0,optTolMinF), 'r-',...
           1:length(fp1), (fp1 - fd1) ./ max(fp1,optTolMinF), 'b-',...
           1:length(fp0), (fp0 - fdBase) ./ max(fp0,optTolMinF), 'g-');
      set(gca,'yscale','log')
      hold on;
      plot([1,max(length(fp0),length(fp1))],[optTol,optTol],'k--');
      hold off
      %}

      figure(1);
      h = plot(1:length(info0.historyFun),(info0.historyFun - fdBase) / max(fdBase,optTolMinF),'b-',...
               1:length(info0.historyDual), (fdBase - info0.historyDual) / max(fdBase,optTolMinF),'r-',...
               1:length(info1.historyDual), (fdBase - info1.historyDual) / max(fdBase,optTolMinF),'b-');
      hold on;
      plot([1,max(length(fp0),length(fp1))],[optTol,optTol],'k--');
      hold off
      set(h(1),'Color',0.7*[1,1,1]);
      set(h,'Linewidth',2);
      set(gca,'Fontsize',16);
      set(gca,'YScale','log')
      xlabel('Iteration');
      ylabel('Relative distance to f^*')

      fprintf('Press <Return> to continue . . . \n')
      pause
   end
end
