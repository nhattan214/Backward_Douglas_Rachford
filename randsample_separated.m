function supp = randsample_separated(N,K,L)
% random sampling K integers from 1--N with spacing at least L 
supp = randi(N-L*(K),K,1);
supp = sort(supp);
supp = supp + (0:K-1)'*L;
end