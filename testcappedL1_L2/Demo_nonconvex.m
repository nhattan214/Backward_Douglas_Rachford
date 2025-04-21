%% Demo for capped L1 - L2 norm


clear; close all
clc

%% parameter settings



matrixtype = 2; 
% 1 for Gaussian
% 2 for DCT
ALL_MEAN_TABLE=[];
for s=1:1

M = 360*s; N = 1280*s;    % matrix dimension M-by-N
K = 40*s;             % sparsity
all_time=[];
all_err=[];
all_iter=[];

for times=1:30
%% construct sensing matrix
switch matrixtype 
    case 1
        A   = randn(M,N); % Gaussian matrix
        A= orth(A')';
    case 2    
        A   = dctmtx(N); % dct matrix
        idx = randperm(N-1);
        A   = A([1 idx(1:M-1)+1],:); % randomly select m rows but always include row 1
        A= orth(A')';
end


%% construct sparse ground-truth 
x_ref = zeros(N,1); % true vector
xs = randn(K,1);
x_ref(randsample(N,K)) = xs;


%% Given lambda, construct b, so that x is a stationary point
lambda = 0.1;

x = x_ref;
[b,y,w,output] = construct_test4L12(A,x,lambda);
% % 
% % 
% % check optimality
norm(lambda * (w - x / norm(x)) + A' * (A * x - b))
max([max(w)-1,min(w) + 1,norm(w(x > 0) - 1),norm(w(x < 0) + 1)])
% 


%% compare L1-L2 solvers
pm.lambda = lambda;
pm.delta = pm.lambda*100;
pm.xg = x; 
pm.reltol  = 1e-6;
pm.maxit = 3000;
pm.height_thres=100;
pmFB = pm; pmFB.delta =  1;

pm_PSAE=pm;
pm_GPPA=pm;


pm_BDR= pm;
pm_BDR.delta=1/20; 
pm_BDR.gamma0=0.15;
pm_BDR.nu=1.4;
pm_BDR.ratio=2;

% initialization
fprintf('Times #%d________________Size #%d_________________________\n', [times,s]);

[xPSAE,outputPSAE] = CS_L1L2_uncon_PSAE_cap(A,b,pm_PSAE,1,50);
[x_BDR,outputBDR] = CS_L1L2_uncon_BDR_cap(A,b,pm_BDR,1);




time=[outputPSAE.time, outputBDR.time];

err=[outputPSAE.err(end),  outputBDR.err(end)];

iter=[ size(outputPSAE.err,2), size(outputBDR.err,2)];


all_time=[all_time; time];
all_err=[all_err; err];
all_iter=[all_iter; iter];
end

MEAN_TIME=mean(all_time);
MEAN_ERR=mean(all_err);
MEAN_ITER=mean(all_iter);
ALL_MEAN=[MEAN_TIME,MEAN_ITER,MEAN_ERR];

ALL_MEAN_TABLE=[ALL_MEAN_TABLE; ALL_MEAN];

end