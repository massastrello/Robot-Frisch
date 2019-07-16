function IDX = CostFcn2(X,n,N,nq,nqd,nqdd)
lambda1 = 1;
lambda2 = 1;
%
Q = reshape(X,[N,3*n]);
%

Wn = computeRegression(Q(:,1:n),Q(:,n+1:2*n),Q(:,2*n+1:end),n,N);
%
%k = 10;
%S = zeros(size(Wn,2),size(Wn,2));
%for i = 1:k
    %nq = randn(N,n);
    %nqd = randn(N,n);
    %nqdd = randn(N,n);
    Wns = computeRegression(Q(:,1:n)+nq,Q(:,n+1:2*n)+nqd,Q(:,2*n+1:end)+nqdd,n,N);
    Ws = Wns-Wn;
    S = ((Ws')*(Ws))/(n*N);
%end
%S = S./k;

%IDX = max(max(S-diag(S))); 
IDX = norm(S-diag(diag(S)),'inf'); %+ norm(diag(S));
%max(max((S-diag(S))));%norm(sort(eig(S))-sort(diag(S)));% max(max((S-diag(S))))/norm(diag(S)) + norm(diag(S));
%figure(3)
%imagesc(S)
%
IDX = IDX + lambda1*cond(Wn) + lambda2*(1/min(svd(Wn)));
%
end

