clc
clear
close all



%% matrixtype choice
%matrix type 1 for Gaussian
%matrix type 2 for DCT


ALL_MEAN_TABLE=[];
for s=1:10

matrixtype = 1; 
M = 360*s; N = 1280*s % matrix dimension M-by-N
K =40*s;             % sparsity


all_time=[];
all_err=[];
all_iter=[];



for times=1:30




%% construct sensing matrix
switch matrixtype
    
    case 1 %Gaussian with noise
        lambda =0.1;
        A = randn(M,N);
        A= orth(A')'; 
        
    case 2 % DCT with noise
        lambda =0.1;
        A   = dctmtx(N); % partial DCT matrix
        idx = randperm(N-1);
        A   = A([1 idx(1:M-1)+1],:); % randomly select m rows but always include row 1
        A= orth(A')';
   
end
    
    
%% construct sparse ground-truth
x_ref = zeros(N,1); % true vector
xs = randn(K,1);
x_ref(randi(N,K,1)) = xs;
x = x_ref;
sigma = 1e-3;
b     = A*x + sigma*randn(M,1);


%% compare L1-L2 solvers
pm.lambda = lambda;
pm.delta = pm.lambda*100;
pm.xg = x_ref;
pm.maxit=2000;
pm.reltol = 1e-6;

pm_BDR=pm;
pm_BDR.delta=1/20; 
pm_BDR.gamma0=0.15;
pm_BDR.nu=1.4;
pm_BDR.ratio=10;


% initialization



fprintf('Times #%d________________Size #%d_________________________\n', [times,s]);


fprintf('Running BDR ...\n');
[x_BDR,outputBDR] = CS_L1L2_uncon_BDR(A,b,pm_BDR,0);


time=[outputBDR.time];

err=[outputBDR.err(end)];

iter=[size(outputBDR.err,2)];


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


