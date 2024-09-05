% Copyright contributors to the Hybrid-BPDN-solver project

function plot_Pareto()

% Determine the filename for caching
prefix = mfilepath(mfilename('fullpath'));
prefix_cache = [prefix,'cache',filesep,'pareto',filesep];
if (~exist(prefix_cache, 'dir'))
   [status,msg,msgID] = mkdir(prefix_cache);
end

filename = [prefix_cache,'figPareto.mat'];
if (exist(filename,'file'))
   % Load the results
   data = load(filename);
   tau0       = data.tau0;
   bNorm      = data.bNorm;
   tauCurve   = data.tauCurve;
   rNormCurve = data.rNormCurve;
   gNormCurve = data.gNormCurve;
   tauRoot    = data.tauRoot;
   rNormRoot  = data.rNormRoot;
   gNormRoot  = data.gNormRoot;

   sigma = 0.1 * bNorm;
else
   % Generate a random problem
   s = RandStream('mt19937ar','Seed',0);
   A = randn(s,1024,2048);
   A = A * diag(1./sqrt(sum(A.^2,1)'));

   k = 450;
   x0 = zeros(2048,1);
   x0(1:k) = randn(s,k,1);
   b = A*x0;
   %v = randn(1024,1); v = v / norm(v,2);
   bNorm = norm(b,2);

   % Find tau_0
   options = struct();
   options.qnType     = 1;
   options.optTol     = 1e-7;
   options.verbosity  = 0;
   options.optTolMinF = 1e-5;
   options.lsType     = 1;
   options.iterations = 5000;
   [x,r,g,info] = solve_rootfinding_v04(A,b,0,options);
   tau0 = norm(x,1);

   % Compute the entire trajectory
   tauCurve = linspace(0,1,101) * tau0;
   rNormCurve = zeros(size(tauCurve));
   gNormCurve = zeros(size(tauCurve));
   x     = zeros(2048,1);
   for i=1:length(tauCurve)
      fprintf('Working on %d/101 . . .\n', i)
      [x,r,g,info] = solver_v04(A,b,0,0,tauCurve(i),[],x,options);

      rNormCurve(i) = sqrt(r'*r);
      gNormCurve(i) = norm(g,Inf);
   end

   % Root finding
   tauRoot   = zeros(1,6);
   rNormRoot = zeros(1,6);
   gNormRoot = zeros(1,6);

   sigma = 0.1 * bNorm;
   tau = 0;
   x = zeros(2048,1);
   for i = 1:6
      [x,r,g,info] = solver_v04(A,b,0,0,tau,[],x,options);
      tauRoot(i) = tau;
      rNormRoot(i) = sqrt(r'*r);
      gNormRoot(i) = norm(g,Inf);

      tau = tau + rNormRoot(i) * (rNormRoot(i) - sigma) / gNormRoot(i);
   end

   % Save the results
   save(filename, 'tau0','bNorm','tauCurve','rNormCurve','gNormCurve', ...
                  'tauRoot','rNormRoot','gNormRoot');
end

% Plotting
figure(1); clf; gca; hold on
tauMax = 1.2 * tau0;
sigmaMax = 1.1*bNorm;

hCurve = plot([tauCurve,tau0,tauMax],[rNormCurve,0,0],'b-');
hSigma = plot([0,tauMax],[sigma,sigma],'r--');
hObj   = plot(tauRoot,rNormRoot,'ro');

hSearch = zeros(1,2);
hTau    = zeros(1,2);
for i=1:3
   hSearch(i) = plot([tauRoot(i),tauRoot(i+1)],[rNormRoot(i),sigma],'r--');
   hTau(i)    = plot([tauRoot(i+1),tauRoot(i+1)],[sigma,rNormRoot(i+1)],'r--');
end
hold off

xlim([0,tauMax]); ylim([0,sigmaMax]);
set(gca,'Fontsize',12);
set(gca,'XTick',[0,tauRoot(end),tau0])
set(gca,'YTick',[0,sigma,bNorm]);
set(gca,'XTickLabel',['0',' ',' ']');
set(gca,'YTickLabel',['0',' ',' ']');
set(hObj,'MarkerSize',3,'LineWidth',3);
set(hCurve,'Linewidth',2);
set(gca,'PlotBoxAspectRatio',	[2.2,1,1]); % 2.4
box on;
