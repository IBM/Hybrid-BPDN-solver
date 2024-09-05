% Copyright contributors to the Hybrid-BPDN-solver project

function table_sparco()

problemIdx = [701, 702, 601, 903, 603, 401, 402];
scaling    = [100, 100,   1, 100, 100,   1,   1];
scales     = [1e-2,1e-3]; % Multiplied by ||b||_2 and scaling to get sigma
legends    = {'Runtime (s)','Rel.gap','Outer iterations'};

xMin        = zeros(length(scales), length(problemIdx));
bNorm       = zeros(1,length(problemIdx));
problemM    = zeros(1,length(problemIdx));
problemN    = zeros(1,length(problemIdx));
problemName = cell(1,length(problemIdx));

for scaleIdx = 1:length(scales)

   % Output table with sparco problems
   fp = 1;
   scale = scales(scaleIdx);

   % ------------------------------------------------------------------
   % Table header
   % ------------------------------------------------------------------
   fprintf(fp,'\\begin{tabular}{lrrcrrcrrcrrcrr}\n');
   fprintf(fp,'&\\multicolumn{2}{c}{SPGL1}&&');
   fprintf(fp,'Orig&Hybr&&Orig&Hybr&&Orig&Hybr&&Orig&Hybr\\\\\n');
   fprintf(fp,'&$10^{-6}$ & $10^{-9}$ &&');
   fprintf(fp,'\\multicolumn{2}{c}{Tol = $10^{-1}$}&&');
   fprintf(fp,'\\multicolumn{2}{c}{Tol = $10^{-2}$}&&');
   fprintf(fp,'\\multicolumn{2}{c}{Tol = $10^{-3}$}&&');
   fprintf(fp,'\\multicolumn{2}{c}{Tol = $10^{-4}$}\\\\\n');
   fprintf(fp,'\\cline{2-3}\\cline{5-6}\\cline{8-9}\\cline{11-12}\\cline{14-15}\n');


   % ------------------------------------------------------------------
   % Compute or retrieve the results
   % ------------------------------------------------------------------
   for i = 1:length(problemIdx)

      table = [];
      % SPGL1
      for methodIdx = [1,2,6]
         %fprintf('Working on problem %d/%.1e, method %d . . .\n', problemIdx(i), scale, methodIdx);
         d = run_sparco(problemIdx(i), scale, methodIdx, 1e-6);
         table(:,end+1) = [d.runtime; d.xNorm1; d.obj; d.rGap; d.info.nNewton];
      end

      su = [];
      for optTol = [1e-1,1e-2,1e-3,1e-4]

         % Hybrid
         runtimePrev = 0; runtime = 0;
         for methodIdx = 3:4
            %fprintf('Working on problem %d/%.1e, method %d (%.1e) . . .\n', problemIdx(i), scale, methodIdx, optTol);
            d = run_sparco(problemIdx(i), scale, methodIdx, optTol);

            runtimePrev = runtime;
            runtime = sum(d.info.runtime);
            table(:,end+1) = [sum(d.info.runtime); d.xNorm1; d.info.info{end}.f; d.info.info{end}.rGap; length(d.info.info)];

            if (d.info.info{end}.rGap > optTol)
               fprintf('Problem %d, method %d, %9.3e vs %9.3e\n', problemIdx(i), methodIdx, d.info.info{end}.rGap, optTol);
            end
         end

         su(end+1) = 100 * (runtimePrev - runtime) / runtimePrev;
      end

      % Get problem information
      try
         p = generateProblem(problemIdx(i));
         problemName{i} = p.info.name;
         problemM(i) = p.sizeA(1);
         problemN(i) = p.sizeA(2);
         bNorm(i) = norm(p.b,2);
      catch exception
         % Hard-coded data to enable table generation without Sparco
         switch problemIdx(i)
             case 401
                problemName{i} = 'srcsep1';
                problemM(i)    = 29166;
                problemN(i)    = 57344;
                bNorm(i)       = 2.2e1;
             case 402
                problemName{i} = 'srcsep1';
                problemM(i)    = 29166;
                problemN(i)    = 86016;
                bNorm(i)       = 2.3e1;
             case 601
                problemName{i} = 'soccer1';
                problemM(i)    = 3200;
                problemN(i)    = 4096;
                bNorm(i)       = 5.5e4;
             case 603
                problemName{i} = 'yinyang';
                problemM(i)    = 1024;
                problemN(i)    = 4096;
                bNorm(i)       = 2.5e1;
             case 701
                problemName{i} = 'blurrycam';
                problemM(i)    = 65536;
                problemN(i)    = 65536;
                bNorm(i)       = 1.3e2;
             case 702
                problemName{i} = 'blurryspike';
                problemM(i)    = 16384;
                problemN(i)    = 16384;
                bNorm(i)       = 2.2e0;
             case 903
                problemName{i} = 'spiketrn';
                problemM(i)    = 1024;
                problemN(i)    = 1024;
                bNorm(i)       = 5.7e1;
         end
      end

      sigma = bNorm(i) * scaling(i) * scale;

      table(3,:) = sqrt(2*table(3,:)); % rNorm
      feasible = (abs(table(3,:) - sigma) / max(1e-3,sigma)) < 1e-5;

      %xMin(scaleIdx,i) = min(table(2,feasible)); % Minimum feasible xNorm
      %table(2,:) = 1e6 * (table(2,:) - xMin(scaleIdx,i)) / xMin(scaleIdx,i);

      table(3,:) = []; % Delete the row: [runtime, xnorm, rgap, nnewton]
      table = table([1,3,4],:); % [runtime, rgap, nnewton]

      % Delete the first column (default spgl1 fails to reach criterion)
      table(:,1) = [];
      feasible(1) = [];

      %fprintf('Problem %d/%e\n', problemIdx(i), scale);
      width = 10;
      fmt = {@formatRuntime, ...
             @formatRelative, ...
             @(v,width) sprintf('%*d',width,v), ...
             @(v,width) sprintf('%*d',width,v)};
      for ii=1:size(table,1)
         fprintf(fp,'%s & ', legends{ii});

         for jj=1:size(table,2)
            % Column separators
            if (jj > 1)
               fprintf(fp, ' & ');
            end
            if ((jj==3) || (jj==5) || (jj==7) || (jj==9))
               fprintf(fp, ' & ');
            end

            formatFun = fmt{ii};
            if (~feasible(jj))
               str = formatFun(table(ii,jj),width);
               fprintf(fp,'\\color{gray}{%s}', str);
            else
               str = formatFun(table(ii,jj),width);
               fprintf(fp,'%s', str);
            end
         end
         fprintf(fp,'\\\\\n');

      end

      %fprintf(fp,'\\hline');
      %fprintf(fp,'\\multicolumn{3}{c}{\\color{blue}{$[$Problem %d$]$}}',problemIdx(i));
      fprintf(fp,'\\multicolumn{3}{l}{\\color{blue}{\\raisebox{1pt}{\\footnotesize$\\blacktriangleright$} Problem %d -- %s}}',problemIdx(i),problemName{i});

      % Speed up
      for jj=1:length(su)
         fprintf(fp,'&&&{\\color{blue}{%.1f}}',su(jj));
      end
      fprintf(fp,'\\\\\n');
      %fprintf(fp,'\\hline\n');

      % End of block
      if (i ~= length(problemIdx))
         fprintf(fp,'\\\\[-7pt]\n');
      end

   end % Problem index

   % ------------------------------------------------------------------
   % Table footer
   % ------------------------------------------------------------------
   fprintf(fp,'\\end{tabular}');

   % Close the function handle
   if (fp ~= 1)
      fclose(fp);
   end

end % Scale


% ------------------------------------------------------------------
% Table with problem information
% ------------------------------------------------------------------
fp = 1;

fprintf(fp,'\\begin{tabular}{lrrrrrrr}\n');
fprintf(fp,'\\hline\n');
fprintf(fp,'Problem & ID & $m$ & $n$ & $\\norm{b}_2$ & scale & $\\norm{x^*_{1e-2}}_1$  & $\\norm{x^*_{1e-3}}_1$\\\\\n');
fprintf(fp,'\\hline\n');
for i=1:length(problemIdx)
   fprintf(fp,'%s & %d & %d & %d & %s & %d & %s & %s\\\\\n',...
           problemName{i},problemIdx(i),problemM(i),problemN(i), ...
           formatNorm(bNorm(i)),scaling(i),formatNorm(xMin(1,i)),formatNorm(xMin(2,i)));
end
fprintf(fp,'\\hline\n');
fprintf(fp,'\\end{tabular}');



function str = formatRuntime(v,width)
   if (v < 100)
      str = sprintf('%*.1f',width,v);
   else
      str = sprintf('%*.0f',width,v);
   end
end

function str= formatDifference(v,width)
   if (v >= 100)
      str = sprintf('%.0f',v);
   else
      str = sprintf('%.1f',v);
   end
end

function str = formatRelative(v,width)
   if (v > 1)
      str = sprintf('%.1f',v);
   else
      str = sprintf('%.0e',v);
      str(end-1) = [];
   end
end

function str = formatNorm(v)
   str = sprintf('%.1e',v);
   str(end-1) = [];
end

end % table_sparco
