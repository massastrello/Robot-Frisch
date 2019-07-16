function IDX = CostFcn(X,n,N)
lambda1 = 1;
lambda2 = 1;
%
Q = reshape(X,[N,3*n]);
%
Wn = computeRegression(Q(:,1:n),Q(:,n+1:2*n),Q(:,2*n+1:end),n,N);
%
%Wns = computeRegression(Q(:,1:n)+randn(N,n),Q(:,n+1:2*n)+randn(N,n),Q(:,2*n+1:end)+randn(N,n),n,N);
%Ws = Wns-Wn;
%S = ((Ws')*(Ws))/(n*N);
%IDX = norm((S-diag(S)),1);
%
IDX = lambda1*cond(Wn) + lambda2*(1/min(svd(Wn)));
%
end

