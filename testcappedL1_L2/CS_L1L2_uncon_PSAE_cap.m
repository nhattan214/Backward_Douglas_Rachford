function [x, output] = CS_L1L2_uncon_PSAE_cap(A,b,pm,activate_restart,restart_step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         min_x .5||Ax-b||^2 + lambda(|x|_1- alpha |x|_2)             %%%
%%%                                                                     %%%
%%% Input: dictionary A, data b, parameters set pm                      %%%
%%%        pm.lambda: regularization paramter                           %%%
%%%        pm.delta: penalty parameter for FBS                          %%%
%%%        pm.maxit: max iterations                                     %%%
%%%        pm.reltol: rel tolerance for FBS: default value: 1e-6        %%%
%%%        pm.alpha: alpha in the regularization                        %%%
%%% Output: computed coefficients x                                     %%%
%%%        output.relerr: relative error of yold and y                  %%%
%%%        output.obj: objective function of x_n:                       %%%
%%%        obj(x) = lambda(|x|_1 - alpha |x|_2)+0.5|Ax-b|^2             %%%
%%%        output.res: residual of x_n: norm(Ax-b)/norm(b)              %%%
%%%        output.err: error to the ground-truth: norm(x-xg)/norm(xg)   %%%
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
% parameter for ADMM
if isfield(pm,'delta');
    delta = pm.delta;
else
    delta = 100 * lambda;
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
    reltol  = 1e-8;
end
if isfield(pm,'alpha');
    alpha = pm.alpha;
else
    alpha = 1;
end

if isfield(pm,'height_thres'); 
    height_thres = pm.height_thres;
else 
    height_thres= 1;
end


%% initialize
x0=zeros(N,1)+10^-6;
x 		= x0;
ell=1;


lambda_bar=0.1;
mu_bar=0.01;
% delta=2*ell*norm(A,norm_type)^2*(lambda_bar^2+0.5)+2*max(norm(A,norm_type)^2*ell*lambda_bar,mu_bar)+0.01;
% delta=2*ell*norm(A,norm_type)^2*(lambda_bar^2+0.5)+2*max(norm(A,norm_type)^2*ell*lambda_bar,mu_bar);

tau_n=1/((2*(5*10^-25))+ell*norm(A'*A)*(2*lambda_bar+1)+2*mu_bar);

xold    = x;
xold2 =x;
% mu_n=tau_n*mu_bar;
t       = 1;
told    = 1;


obj         = @(x) .5*norm(A*x-b)^2 + lambda*(norm(x,1)-alpha*norm(x));
output.pm   = pm;




temp=(tau_n).*A'*b;
temp2=(tau_n).*A'*A;



step_size=lambda*tau_n;

total_time=0;

for it = 1:maxit
    tic;

    u_n=xold+((told-1)/t)*lambda_bar.*(xold-xold2);

    v_n=xold+tau_n*mu_bar.*(xold-xold2);

    told = t;
    t = (sqrt(4*t^2+1) + 1)/2;
    if activate_restart==1
        if mod(it,restart_step)==0
            told=1;
            t=1;
        end
    end
   
    x=prox_CapL1(v_n-temp2*u_n+temp+step_size.*(xold/norm(xold,2)),height_thres,step_size);


    time_iter=toc;
    % stop conditions & outputs
 
    relerr=norm(x-xold)/norm(xold);
    xold2=xold;
    xold=x;
    %     obj_current=obj(x);
    %     relerr=abs(obj_prev-obj_current);

    residual    = norm(A*x - b)/norm(b);

    output.relerr(it)   = relerr;
    output.obj(it)      = obj(x);

    
    output.res(it)      = residual;
    output.err(it)      = norm(x - xg)/norm(xg);
    total_time = total_time + time_iter;
    output.time_cum(it)     = total_time;
        
    if relerr < reltol && it>2
        fprintf('Stopped at iter %d by relative error\n',it);
        output.Stopped_at=it;
        break;
        
    end
end

if it==maxit
    fprintf('Stopped by max iteration');
    output.Stopped_at=it;

end
% endtime=toc(start_time)
output.time = total_time;

end




