% Copyright contributors to the Hybrid-BPDN-solver project

function plot_problem_duality()

nInstances     = 10;
sparsityLevels = 7;
tauLevels      = 6;
noiseLevels    = 4;

nJobs = sparsityLevels * tauLevels * noiseLevels;

jobs = 1:nJobs;

compare1 = 1;
compare2 = 2;

colorLevels = tauLevels; % sparsityLevels;

% Selected points for each of the settings
settings = [1, 5, 2,  4; ...  # muIdx=1
            2, 5, 1,  1; ...  # muIdx=2 (1,6,2,3)
            7, 6, 2, 10; ...  # muIdx=3
            5, 6, 3,  7];

for muIdx=1:4 % 1-4
   s = 0;

   figure(1); clf; gca; hold on;
   set(gca,'Fontsize',16);
   v = hsv(colorLevels);
   xlabel('Runtime original dual (s)');
   ylabel('Speed up');

   for i=1:length(jobs)
      k = jobs(i) - 1;

      sparsityIdx = mod(k,sparsityLevels) + 1;
      k = floor(k / sparsityLevels);
      tauIdx = mod(k,tauLevels) + 1;
      k = floor(k / tauLevels);
      noiseIdx = mod(k,noiseLevels) + 1;

      % 1 = original without
      % 2 = original with
      % 3 = hybrid without
      % 4 = hybrid with
      r = run_problem_duality(muIdx,sparsityIdx, tauIdx, noiseIdx, nInstances);

      % Plot speedup versus runtime
      x = r.runtime(:,compare1);
      y = r.runtime(:,compare1) ./ r.runtime(:,compare2);

      idx1 = (r.stat(:,compare1) == 4);
      idx2 = (r.stat(:,compare2) == 4);

      % Plot the results
      idx = (idx1 & idx2);
      h1 = plot(x(idx),y(idx),'b.'); set(h1,'color',0.45*[1,1,1]);
      idx = (~idx1 &  idx2);
      h2 = plot(x(idx),y(idx),'r.');
      idx = (idx1 & ~idx2);
      h3 = plot(x(idx),y(idx),'b.');
      idx = (~idx1 & ~idx2);
      h4 = plot(x(idx),y(idx),'k.');

      if ((sparsityIdx == settings(muIdx,1)) && ...
            (tauIdx      == settings(muIdx,2)) && ...
            (noiseIdx    == settings(muIdx,3)))
         idx = settings(muIdx,4);
         h = plot(x(idx),y(idx),'bo');
         set(h,'markersize',12,'linewidth',1.5);
      end

      % Keep track of the total runtime
      s = s + sum(r.runtime,1);
   end
   disp(s);
   hold off; box on;

   fprintf('Press <Return> to continue . . .\n');
   pause;
end
