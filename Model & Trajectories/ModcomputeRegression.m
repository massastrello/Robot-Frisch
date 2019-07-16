function Wn = ModcomputeRegression(q,dq,ddq)
N = min([size(q,1) size(dq,1) size(ddq,1)]);
Wn = [];
for i = 1:N
    Wi = ModmakeWi(q(i,:),dq(i,:),ddq(i,:));
    Wn = [Wn;Wi];
end
end