[x1,fval1] = particleswarm(@(Xi) CostFcn(Xi,n,N),3*n*N,LB,UB,options);
T1 = reshape(x1,[N,3*n]);
T1 = T1';%[zeros(3*n,1) T' zeros(3*n,1)]; % add a column of zero as initial point of q, dq, ddq
Q1 = T1';
q1 = Q1(:,1:n);
qd1 = Q1(:,n+1:2*n);
qdd1 = Q1(:,2*n+1:3*n);

nt = 1000;
S = zeros(9,9);
S1 = zeros(9,9);
Srand = zeros(9,9);
for i = 1:nt
rng(i);
nqrand = randn(N,n);%nq;%
rng(i + 2*nt)
nqdrand = randn(N,n);%nqd;%
rng(i + 4*nt)
nqddrand = randn(N,n);%nqdd;%

rng(i + 6*nt)
qrand = randn(N,n);
rng(i + 8*nt)
qdrand = randn(N,n);
rng(i + 10*nt)
qddrand = randn(N,n);
%
Wn = computeRegression(q,qd,qdd,n,N);
Wn1 = computeRegression(q1,qd1,qdd1,n,N);

Wnrand = computeRegression(qrand,qdrand,qddrand,n,N);
Wns = computeRegression(q+nqrand,qd+nqdrand,qdd+nqddrand,n,N);
Wns1 = computeRegression(q1+nqrand,qd1+nqdrand,qdd1+nqddrand,n,N);
Wnsrand = computeRegression(qrand+nqrand,qdrand+nqdrand,qddrand+nqddrand,n,N);

Ws = Wns-Wn;
Ws1 = Wns1-Wn1;
Wsrand = Wnsrand-Wnrand;

S  = S + ((Ws')*Ws)/(n*N);
S1  = S1 + ((Ws1')*Ws1)/(n*N);
Srand  = Srand + ((Wsrand')*Wsrand)/(n*N);
end
S = S./nt;
S1 = S1./nt;
Srand = Srand./nt;
%[x,fval] = particleswarm(@(Xi) CostFcn2(Xi,n,N,nq,nqd,nqdd),3*n*N,LB,UB,options);

figure()
subplot(131)
imagesc(S-diag(diag(S)))
subplot(132)
imagesc(S1-diag(diag(S1)))
subplot(133)
imagesc(Srand-diag(diag(Srand)))

E = (abs(S1-diag(diag(S1)))-abs(S-diag(diag(S))));
e_mean = mean(reshape(E,length(E)^2,1))
figure()
imagesc(E)
%colormap('gray')
colorbar
