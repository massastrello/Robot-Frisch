function Wn = computeRegression(q,dq,ddq,n,N)
Wn = zeros(N*n,size(makeWi([0 0],[0 0],[0 0]),2));
j = 1;
for i = 1:N
   Wn(j:n+j-1,:) = makeWi(q(i,:),dq(i,:),ddq(i,:));
   j = j + n;
end
end