% Copyright contributors to the Hybrid-BPDN-solver project

function plot_Lipschitz()


% Determine the filename for caching
prefix = mfilepath(mfilename('fullpath'));
prefix_cache  = [prefix,filesep,'cache',filesep];
if (~exist(prefix_cache, 'dir'))
   [status,msg,msgID] = mkdir(prefix_cache);
end

% Set the global random stream
stream = RandStream.getGlobalStream();
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

% Problem settings
m = 200; n = 500;

A = randn(m,n);
A = A * diag(1./sqrt(sum(A.^2,1)'));

gammaVals = linspace(0,0.1,101);

mu = zeros(1,length(gammaVals));
sg = zeros(1,length(gammaVals));

for i=1:length(gammaVals)
   B = random_walk_matrix(A,1-gammaVals(i));
   G = B'*B; G = G - diag(diag(G));
   S = svd(B);
   mu(i) = max(abs(G(:)));
   sg(i) = S(1);
end
h = plot(gammaVals,sg);
set(h,'Linewidth',2);
set(gca,'Fontsize',22);
xlabel('\gamma'); ylabel('\sigma_1');

fprintf('Press <Return> to continue . . .\n');
pause;

% Mutual coherence and largest singular value
m = 200; n = 500;
filename = [prefix_cache,'FigLipschitz.mat'];
if (exist(filename,'file'))
   data = load(filename);
   mu = data.mu;
   sg = data.sg;
else
   mu = zeros(1,1000);
   sg = zeros(1,1000);
   for i=1:1000
      disp(i);
      A = randn(m,n);
      A = A * diag(1./sqrt(sum(A.^2,1)'));
      G = A'*A; G = G - diag(diag(G));
      S = svd(A);
      mu(i) = max(abs(G(:)));
      sg(i) = S(1);
   end
   save(filename,'mu','sg');
end
h = plot(mu,sg,'b*');
set(gca,'Fontsize',22);
xlabel('Mutual coherence'); ylabel('\sigma_1');

fprintf('Press <Return> to continue . . .\n');
pause;


% Gram matrix
m = 200; n = 2000;
gammaValues = [0.01, 0.005];
for i=1:length(gammaValues)
   gamma = gammaValues(i);
   A = randn(m,n);
   B = random_walk_matrix(A,1-gamma);
   pcolor(B'*B),shading flat
   colormap(jet);
   set(gca,'Fontsize',22); box on;

   if (i ~= length(gammaValues))
      fprintf('Press <Return> to continue . . .\n');
      pause;
   end
end

% Restore the original global random stream
RandStream.setGlobalStream(stream);
