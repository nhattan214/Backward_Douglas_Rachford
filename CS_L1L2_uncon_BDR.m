function [z, output] = CS_L1L2_uncon_BDR(A,b,pm,heuristic_on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         min_x .5||Ax-b||^2 + lambda(|x|_1- alpha |x|_2)             %%%
%%%    phi = 0.5||Ax-b||^2; g = lambda |x|_1; h = -alpha*lambda*|x|_2   %%%
%%%                                                                     %%%
%%% Input: dictionary A, data b, parameters set pm                      %%%
%%%        pm.lambda: regularization paramter                           %%%
%%%        pm.delta: penalty parameter for DRDC                         %%%
%%%        pm.maxit: max iterations                                     %%%
%%%        pm.reltol: rel tolerance for DRDC: default value: 1e-6       %%%
%%%        pm.alpha: alpha in the regularization,default value: 1       %%%
%%% Output: computed coefficients z (shrinkage operator result)         %%%
%%%        output.relerr: relative error of yold and y                  %%%
%%%        output.obj: objective function of x_n:                       %%%
%%%        obj(x) = lambda(|x|_1 - alpha |x|_2)+0.5|Ax-b|^2             %%%
%%%        output.res: residual of x_n: norm(Ax-b)/norm(b)              %%%
%%%        output.err: error to the ground-truth: norm(z-xg)/norm(xg)   %%%
%%%        output.time: computational time                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N]       = size(A); 
% start_time  = tic; 

%% parameters
if isfield(pm,'lambda'); 
    lambda = pm.lambda; 
else
    lambda = 1e-5;  % default value
end

% maximum number of iterations
if isfield(pm,'maxit'); 
    maxit = pm.maxit; 
else 
    maxit = 5*N; % default value
end
% initial guess
if isfield(pm,'x0'); 
    x0 = pm.x0; 
else 
    x0 = zeros(N,1); % initial guess
end
if isfield(pm,'xg'); 
    xg = pm.xg; 
else 
    xg = x0;
end
if isfield(pm,'reltol'); 
    reltol = pm.reltol; 
else 
    reltol  = 1e-6; 
end
if isfield(pm,'alpha'); 
    alpha = pm.alpha;
else 
    alpha = 1;
end

if isfield(pm,'delta'); 
    delta = pm.delta;
else 
    delta = 1;
end

if isfield(pm,'nu'); 
    nu = pm.nu; 
else
    nu = 1;
end

if isfield(pm,'gamma0'); 
    gamma_0 = pm.gamma0;
else 
    gamma_0 = 1;
end

if isfield(pm,'ratio'); 
    ratio = pm.ratio;
else 
    ratio = 1;
end

%% initialize
x 		= x0; 
z 		= x0;
w       = x0;

ell=norm(A*A');

if heuristic_on==1

gamma_0=min(1/ell,(sqrt(8*(2-nu)*ell^2)/(4*ell^2)))-(10^-100)-0.1;
gamma = gamma_0*ratio; 
   
else
gamma=min(1/ell,(sqrt(8*(2-nu)*ell^2)/(4*ell^2)))-(10^-10);
end


obj         = @(x) .5*norm(A*x-b)^2 + lambda*(norm(x,1)-alpha*norm(x));
output.pm   = pm;

% method 1
D_inv = inv(gamma*A'*A + eye(N));
% method 2
L = chol(speye(M) + gamma*(A*A'), 'lower');
L = sparse(L);
U = sparse(L');

Atb   = A'*b;

q_handle=@(q) q - gamma*(A'*(U \ ( L \ (A*q) )));
Ahandle  = @(x) D_inv*x;

total_time = 0;


y=rand(N,1);

for it = 1:maxit
    output.err(it)      = norm(z - xg)/norm(xg);
    
    tic;
    % select k size
    if it==1
        zn=z;
    else
        zn=zn;
    end
    
    %x subproblem
    q = gamma*Atb + w;
    xn = q - gamma*(A'*(U \ ( L \ (A*q) )));
   

   
   


    % y subproblem
    u=y./delta+zn;
    y=min(lambda/norm(u),delta).*u;

    % z subproblem
    zn = shrink(2*xn-w+gamma*y,lambda*gamma);

    % x-update
    wn = w +  nu*(zn - xn);
    time_iter = toc;
   
  
    % stop conditions & outputs

    relerr      = norm(xn - x)/norm(x);
    residual    = norm(A*zn - b)/norm(b);

    if heuristic_on ==1
    if (norm(xn - x) >= 1e3/it || norm(x) > 1e10) && gamma > gamma_0 % decrease gamma to guarantee convergence
       gamma = max(gamma/2,gamma_0*0.9999);
    end
    else
       gamma = gamma;
    end
    
    x = xn;
    w = wn;
    z = zn;
    
    
    output.relerr(it)   = relerr;
    output.obj(it)      = obj(z);
    output.res(it)      = residual;
    total_time = total_time + time_iter;
    output.time_cum(it)     = total_time;
    if relerr < reltol && it > 2  
        break;
    end
end
output.time = total_time;
end


function z = shrink(x, r)
    z = sign(x).*max(abs(x)-r,0);
end

