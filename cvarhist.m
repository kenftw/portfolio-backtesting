% CVaR calculation for scenario returns
% w is in (n x 1) vector and r is a (T x n) vector, where
% n is the number of assets and T the number of simulations

function y = cvarhist(w,r,alpha) 
[T,n] = size(r);
rp = sort(r*w); % sort portfolio returns in ascending order (e.g.,1,2,3)
VaR = floor(alpha*T); % VaR
y = mean(rp(1:VaR,1));
