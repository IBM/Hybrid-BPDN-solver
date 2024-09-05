% Copyright contributors to the Hybrid-BPDN-solver project

function [A,b,tau] = generate_problem_duality(idxSparsity, idxTau, idxNoise, index)

sparsityLevels = [50,100,150,200,250,300,350];
tauMultipliers = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2];
noiseLevels    = [0, 0.0001, 0.001, 0.01];

s = RandStream('mt19937ar','Seed',index + 13 * idxSparsity + 71 * idxTau + 191 * idxNoise);

% Sensing matrix
m = 1000; n = 2000;
A = randn(s,m,n) / sqrt(m);

% Create the sparse signal
k = sparsityLevels(idxSparsity);
x0 = zeros(n,1);
x0(1:k) = randn(s,k,1);

% Observation
b = A*x0;
s = randn(s,m,1);
b = b + (noiseLevels(idxNoise) * norm(b,2) / norm(s,2)) * s;

tau = tauMultipliers(idxTau)*norm(x0,1);
