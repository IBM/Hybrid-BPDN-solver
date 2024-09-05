% Copyright contributors to the Hybrid-BPDN-solver project

function plot_mutual_coherence()

m = 64; n = 2048;

s = RandStream('mt19937ar','Seed',0);
A = randn(s,m,n);
A = A * spdiags(1./sqrt(sum(A.^2,1))',0,n,n);
v = randn(s,m,1); v = v / norm(v,2);

x = linspace(0,1,201);
x = (x(1:end-1) + x(2:end)) / 2;

colors = [0,0.8,0; ...
          0,0,1; ...
          1,0,0];

count = zeros(3,length(x));

for i=1:3
   switch i
      case 1
         M = A;
      case 2
         M = random_walk_matrix(A,1-0.005); % 6
      case 3
         M = random_walk_matrix(A,1-0.002); % 2
   end

   G = M'*M; G = G + diag(inf(n,1));
   G = sort(abs(G(:)));
   G = G(1:end-n);

   c= hist(G,x);
   c = c / (sum(c) * (x(2) - x(1)));

   c(2:end-1) = 0.2 * c(1:end-2) + 0.6 * c(2:end-1) + 0.2 * c(3:end);

   count(i,:) = c;
   G = [];
end

idx = [1,2,3];
h = plot(x,count(idx,:));
for j=1:length(idx)
   i = idx(j);
   set(h(j),'Color',colors(i,:));
end
set(h,'LineWidth',1.2);
set(gca,'PlotBoxAspectRatio',	[3.25,1,1]);
set(gca,'Fontsize',8);
box on;

legend('Normal distribution','\gamma = 0.005','\gamma = 0.002');
xlabel('Pairwise mutual coherence'); ylabel('Density');
