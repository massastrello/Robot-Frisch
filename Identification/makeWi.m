%% Last Edit 07/05/2017

function Wi = makeWi(q,dq,ddq)
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
Y;
%
if isequal(AI,[])
    Wi = W;
    Wi(:,s) = [];
else
    Wii = W;
    Wii(:,s) = [];
    Wc;
    Wi = [Wii(:,AI),Wci];
end