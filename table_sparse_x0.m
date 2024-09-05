% Copyright contributors to the Hybrid-BPDN-solver project

function table_sparse_x0()

nTrials = 50;
sparsityLevels = 50:25:400;

i1   = zeros(3,length(sparsityLevels));
i2   = zeros(3,length(sparsityLevels));
t1   = zeros(3,length(sparsityLevels));
t2   = zeros(3,length(sparsityLevels));
s1   = zeros(3,length(sparsityLevels));
s2   = zeros(3,length(sparsityLevels));
f1   = zeros(3,length(sparsityLevels));
f2   = zeros(3,length(sparsityLevels));
gap1 = zeros(3,length(sparsityLevels), nTrials);
gap2 = zeros(3,length(sparsityLevels), nTrials);

for idx = 1:length(sparsityLevels)
   sparsity = sparsityLevels(idx);
   for type=1:3
      v1 = run_sparse3_x0(1024,2048,sparsity,type,1,nTrials);
      v2 = run_sparse3_x0(1024,2048,sparsity,type,2,nTrials);

      % Iterations
      i1(type,idx) = mean(v1(:,1));
      i2(type,idx) = mean(v2(:,1));

      % Runtime
      t1(type,idx) = mean(v1(:,2));
      t2(type,idx) = mean(v2(:,2));

      % Successful solves
      s1(type,idx) = sum(v1(:,3) == 4);
      s2(type,idx) = sum(v2(:,3) == 4);

      % Objective value
      f1(type,idx) = min(v1(:,4));
      f2(type,idx) = min(v2(:,4));

      % Duality gap
      gap1(type,idx,:) = (v1(:,4) - v1(:,5)) ./ max(v1(:,4), 1e-3);
      gap2(type,idx,:) = (v2(:,4) - v2(:,5)) ./ max(v2(:,4), 1e-3);

      % Trajectory line search
      %v3 = run_sparse2_x0(1024,2048,sparsity,type,3,nTrials);
      %s3(type,idx) = sum(v3(:,3) == 4);
   end
end


% =====================================================================
% Table 1 -- Recovery rate, runtime, and iterations
% =====================================================================
fp = 1;

% ---------------------------------------------------------------------
% Table header
% ---------------------------------------------------------------------
fprintf(fp,'\\begin{tabular}{rcrrrcrrrcrrr}\n');
fprintf(fp,['$k$ && ',...
            '\\multicolumn{3}{c}{Runtime original (s)}&&',...
            '\\multicolumn{3}{c}{Runtime hybrid (s)}&&',...
            '\\multicolumn{3}{c}{Speed up (\\%%)}\\\\\n']);
fprintf(fp,'\\cline{1-1}\\cline{3-5}\\cline{7-9}\\cline{11-13}\\cline{15-17}\n');


% ---------------------------------------------------------------------
% Table body
% ---------------------------------------------------------------------
for i=1:length(sparsityLevels)
   fprintf(fp,'%3d', sparsityLevels(i));
   fprintf(fp,'&');
   for j=1:3, fprintf(fp,'& %5.2f', t1(j,i)); end
   fprintf(fp,'&');
   for j=1:3, fprintf(fp,'& %5.2f', t2(j,i)); end
   fprintf(fp,'&');
   for j=1:3, fprintf(fp,'& %1.f', 100*(1-t2(j,i) / t1(j,i))); end
   fprintf(fp,'\\\\\n');
end

% ---------------------------------------------------------------------
% Table footer
% ---------------------------------------------------------------------
fprintf(fp,'\\end{tabular}');
fprintf(fp,'\n');


% =====================================================================
% Table 2 -- Success rate
% =====================================================================
fp = 1;

fprintf(fp,'\\begin{tabular}{rcrrrcrrr}\n');

for i=1:length(sparsityLevels)
   fprintf(fp,'%3d', sparsityLevels(i));
   fprintf(fp,'&');
   for j=1:3, fprintf(fp,'& %3d', s1(j,i)); end
   fprintf(fp,'&');
   for j=1:3, fprintf(fp,'& %3d', s2(j,i)); end
   %fprintf(fp,'&');
   %for j=1:3, fprintf(fp,'& %3d', s3(j,i)); end
   fprintf(fp,'\\\\\n');
end

fprintf(fp,'\\end{tabular}');
fprintf(fp,'\n');



% =====================================================================
% Table 3 -- Objective value
% =====================================================================
fp = 1;

fprintf(fp,'\\begin{tabular}{rcrrrcrrr}\n');

for i=1:length(sparsityLevels)
   fprintf(fp,'%3d', sparsityLevels(i));
   fprintf(fp,'&');
   for j=1:3, fprintf(fp,'& %3d', f1(j,i)); end
   fprintf(fp,'&');
   for j=1:3, fprintf(fp,'& %3d', f2(j,i)); end
   fprintf(fp,'\\\\\n');
end

fprintf(fp,'\\end{tabular}');
fprintf(fp,'\n');



% =====================================================================
% Table 4 -- Recovery rate, runtime, iterations, and success rate
% =====================================================================
fp = 1;

% ---------------------------------------------------------------------
% Table header
% ---------------------------------------------------------------------
fprintf(fp,'\\begin{tabular}{rcrrrrrcrrrrrcr}\n');
fprintf(fp,['$k$ && ',...
            '\\multicolumn{3}{c}{\\ \\ Runtime original (s)}& rel.gap&\\ding{51}&&',...
            '\\multicolumn{3}{c}{\\ \\ Runtime hybrid (s)}&rel.gap&\\ding{51}&&',...
            '\\multicolumn{1}{c}{(\\%%)}\\\\\n']);
fprintf(fp,'\\cline{1-1}\\cline{3-7}\\cline{9-13}\\cline{15-15}\n');


% ---------------------------------------------------------------------
% Table body
% ---------------------------------------------------------------------
for i=1:length(sparsityLevels)
   % Sparsity level
   fprintf(fp,'%3d', sparsityLevels(i));

   % ====== Original method ======
   fprintf(fp,'&');

   % Run times
   for j=1:3, fprintf(fp,'& %5.2f', t1(j,i)); end

   % Relative duality gap (original)
   x = [squeeze(gap1(1,i,:)); squeeze(gap1(2,i,:)); squeeze(gap1(3,i,:))];
   n = sum(x <= 1e-6);
   x = median(x);
   x = sprintf('%.1e', x); x = [x(1:5),x(7)];
   fprintf(fp,'& %s', x);

   % Percentage of successful solves
   fprintf(fp,'&%3d', round((100 / nTrials)*mean(s1(1:3,i))));

   % ====== Hybrid method ======
   fprintf(fp,'&');

   % Run times
   for j=1:3, fprintf(fp,'& %5.2f', t2(j,i)); end

   % Relative duality gap (hybrid)
   x = [squeeze(gap2(1,i,:)); squeeze(gap2(2,i,:)); squeeze(gap2(3,i,:))];
   n = sum(x <= 1e-6);
   x = median(x);
   x = sprintf('%.1e', x); x = [x(1:5),x(7)];
   fprintf(fp,'& %s', x);

   % Percentage of successful solves
   fprintf(fp,'&%3d', round((100 / nTrials)*mean(s2(1:3,i))));

   % ====== Comparison ======
   fprintf(fp,'&');

   % Speed up (average of three speed-up values)
   s = zeros(1,3);
   for j=1:3
      s(j) = 100*(1-t2(j,i) / t1(j,i));
   end
   fprintf(fp,'& %1.f', mean(s));



   fprintf(fp,'\\\\\n');
end

% ---------------------------------------------------------------------
% Table footer
% ---------------------------------------------------------------------
fprintf(fp,'\\end{tabular}');
fprintf(fp,'\n');
