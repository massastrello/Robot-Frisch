function Wi = ModmakeWi(q,dq,ddq)
global s
global AI

g = 9.81;

q1 = q(1);
q2 = q(2);

qd1 = dq(1);
qd2 = dq(2);

qdd1 = ddq(1);
qdd2 = ddq(2);

%
Y
Wi = W;
Wi(:,s) = [];

% If the absulutely identifi  
if isequal(AI,[])
    % Delete Null Columns
else
Wi(:,AI) = [];
end
   