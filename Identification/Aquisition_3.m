%{
q = qint';
qd = qdint';
qdd = qddint';
%}
N = length(q);
%% ACQUISITION OF DATA AND FRISCH
bn = .001;
tn = 1;
scalingF = 10;
hh=rand;
%% Simulate the behaviour of the robot
par = size(gammaR,2); 
% Add some (bounded) random noise to [q] to simulate measurment errors
% introduced by the sampling process and reconstruction of ddq and dq from
qS = zeros(N,n);
qdS = qS;
qddS = qdS;

coeff_ref2;
%% Select one of the two ways to compute the regression matrix
rng(1+hh)
nq = bn*randn(N,n);
rng(2+hh)
ndq = 2*bn*randn(N,n);
rng(3+hh)
nddq = 2*bn*randn(N,n);

qS = q + nq;
qdS = qd + ndq;
qddS = qdd + nddq;

Wnp = computeRegression(q,qd,qdd,n,N);
Wns = computeRegression(qS,qdS,qddS,n,N); 
%% Declare the parameters vector
par = size(Wnp,2); % Number of parameters to be estimated
%Inertial Parameters Vector
%base parameters
%pi = zeros(par,1);
%for i = 1:par
%    pi(i) = i;
%end
%% Check if, with nominal measurement points and torque, the simplex degenerates in one point
%
disp('* CHECK: Test of Estimation Procedure ...');
We_c = [Wnp, Wnp*gammaR_ref]; %add the torques to the regression matrix in order to perform the OLS
Sigma_c = (We_c'*We_c)/N;
alpha_c = OLS(Sigma_c,par+1,-1);
%
%%
Check = 1;
while (i <= size(alpha_c,2)) && Check
    x = norm([gammaR_ref;-1]-alpha_c(:,i));
    if x<1e-4
        i = i + 1;
    else
        Check = 0;
    end
end
if Check
    msg = '* The OLS found the exact solution with nominal data: the procedure is working!';
    msgg = '   (The simplex degenerates to a single point)';
else
    msg = '* Something went wrong while trying to perform the OLS: the solution is not exact with nominal data';
    msgg = '';
    return;
end
disp(msg);
disp(msgg);
disp('');
%}
%%%%%%%%%%%%%%%%%%%% Perform the OLS to find the vertexes of the symplex in the parameter's space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rng(4+hh)
ntau = tn*randn(size(Wns,1),1);
tau =  Wnp*gammaR_ref + ntau;
tau = tau*scalingF;
We = [Wns, tau]; %add the torques to the regression matrix in order to perform the OLS
Sigma = (We'*We)/N;
alpha = OLS(Sigma,par+1,-1);
alpha = alpha/scalingF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% each row of theta contains the bound in which one of the estimated parameter
% lives and its true value
theta = zeros(par,3);
for i = 1:par
    theta(i,:) = [min(alpha(i,:)) max(alpha(i,:)) gammaR_ref(i)];
end
%% Simplex Check
% Check if the real parameters are inside the simplex optained performing the OLS 
% with the experimentally aquired data
disp('* CHECK: Correctness of the Simplex...');
Check = 1;
Count = 0;
min_e = zeros(1,par);
err = [];
for i = 1:par
    if (theta(i,3) <= theta(i,2)) && (theta(i,3) >= theta(i,1))
        obj_f = [];
        for j = 1:par+1
            obj_f = [obj_f,abs(alpha(i,j)-theta(i,3))];
        end
        min_p(i) = min(obj_f);
    else
        Check = 0;
        if theta(i,3) > theta(i,2)
            min_e(i) = -abs(theta(i,3) - theta(i,2));
        else
            min_e(i) = -abs(theta(i,1) - theta(i,3));
        end
        warning('* pi(%d) is not inside the siplex\nerr = %f\n',i,min_e(i));
        err = [err,min_e(i)];
        Count = Count + 1;
    end
end

if Check == 1
    disp('* All the parameter are inside the simplex'); 

else
    disp(['* ',num2str(Count),' parameter are outside the simplex (Mean error: ',num2str(mean(err)),')']);
end
LSvsFR = zeros(par,3);
for i = 1:par 
    LSvsFR(i,1) = i;
    LSvsFR(i,3) = min(abs(gammaR_ref(i)*ones(1,par) - alpha(i,1:par))) - (abs(gammaR_ref(i) - alpha(i,par+1)));
    if min(abs(gammaR_ref(i)*ones(1,par) - alpha(i,1:par))) - (abs(gammaR_ref(i) - alpha(i,par+1))) < 0
        LSvsFR(i,2) = 1;
    end
end
%disp('Comparison between classical LS Technique and the Frish Scheme:');
%disp('each element of the second column is 1 if there exists a point belonging to the simplex in the parameters space,'); 
%disp('found under conditions of the Frisch Scheme which is more close to the real value of  meter than the least square');
x = sum(LSvsFR(:,2));
disp(['For ',num2str(x),'/',num2str(par),' parameters the Frisch Scheme "got" a closer solution (vertex) than LS (see the variable ''LSvsFR'')']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��%%%%%%%%%%%
%% Add to theta 2 columns with mean and median of simplex
Mean = zeros(par,1);
Median = Mean;
for i = 1:par
    Mean(i,1) = mean(alpha(i,:));
    Median(i,1) = median(alpha(i,:));
end
theta = [theta,Mean,Median];
[gammaR_ref,min(alpha(1:end-1,:)')',max(alpha(1:end-1,:)')',alpha(1:end-1,10)]
clear Check msg msg1