function [Array] = get_range(N,K)

M = rem(N,K);
Array = 1:K:N+K-M;
Array(2:end) = Array(2:end) + N - Array(end);
