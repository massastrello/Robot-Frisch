%% ACQUISITION OF DATA AND FRISCH
%{
bn = 0.05;
tn = 0.02;
scalingF = 10;
%% Simulate the behaviour of the robot
par = size(gammaR,2); 
% Add some (bounded) random noise to [q] to simulate measurment errors
disp('* Sampling Data...');
n_int = length(qint);
qint_noisy = qint - bn*ones(n,n_int) + 2*bn*rand(n,n_int); 
disp('  Done.');
%
figure()
plot(T,qint_noisy,'LineWidth',2)
hold on
plot(T,qint,'k','LineWidth',2)
hold off
% differentiate
qdint_noisy = robustDiff(qint_noisy',ts,101);
figure()
plot(T,qdint_noisy,'LineWidth',2)
hold on
plot(T,qdint,'k','LineWidth',2)
hold off
%}
% introduced by the sampling process and reconstruction of ddq and dq from
qS = zeros(N,n);
dqS = qS;
ddqS = dqS;

coeff_ref;
%% Select one of the two ways to compute the regression matrix
% build up the regression matrix with the reference measurments points or with the noisy sampling
Wnp = computeRegression(q,dq,ddq,n,N);
%Wns = computeRegression(qS,dqS,ddqS,n,N); 
%% Declare the parameters vector
p = size(Wnp,2); % Number of parameters to be estimated

%% Check if, with nominal measurement points and torque, the simplex degenerates in one point
%
disp('* CHECK: Test of Estimation Procedure ...');
We_c = [Wnp, Wnp*gammaR_ref]; %add the torques to the regression matrix in order to perform the OLS
Sigma_c = (We_c'*We_c)/N;
alpha_c = OLS(Sigma_c,p+1,-1);
%
%%
Check = 1;
while (i <= size(alpha_c,2)) && Check
    x = norm([pi;-1]-alpha_c(:,i));
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
tau =  -tn*ones(size(Wns,1),1) + 2*tn*randn(size(Wns,1),1) + Wns*gammaR_ref;
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

clear Check msg msg1