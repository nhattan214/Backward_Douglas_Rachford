function p = prox_CapL1(x,gamma,lambda)
%function p = prox_CapL1(x,gamma,theta)
%
% This procedure computes the proximity operator of the  capped-l1 (or
% truncated L1) penalty [1].
%
% f(x) =  gamma*min(theta*|x|, 1).  
% gamma is the height threshold
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x      - ND array
%  gamma  - positive, scalar or ND array with the same size as 'x'
%  theta  - positive, scalar or ND array with the same size as 'x'
%
%  REFERENCE
% ===========
% [1] T. Zhang. Multi-stage Convex Relaxation for Learning with Sparse Regularization.
%     Neural Information Processing Systens, 1929-1936, 2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (21-10-2019)
% Author  : Emmanuel Soubies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x''')
end
theta=1/gamma;
if any( theta(:) <= 0 ) || ~isscalar(theta) && any(size(theta) ~= size(x))
    error('''theta'' must be positive and either scalar or the same size as ''x''')
end
%-----%

Gamma=gamma*lambda;

p= sign(x).*(abs(x)-Gamma.*theta).*(abs(x) > Gamma.*theta) .* (abs(x)<(1./theta +Gamma.*theta/2)).* (Gamma.*theta.^2 <2) ...
   + x .* (Gamma.*theta.^2 <2).* (abs(x)>=(1./theta +Gamma.*theta/2)) ...
   + x .* (abs(x)>sqrt(2*Gamma)) .* (Gamma.*theta.^2 >= 2);                                                                       
end

