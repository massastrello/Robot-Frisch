%
q = qint';
qd = qdint';
qdd = qddint';
%}
N = length(q);
%% ACQUISITION OF DATA AND FRISCH
bn = .001;
tn = 1;
scalingF = 1;
hh=10;
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
%
%% Check if, with nominal measurement points and torque, the simplex degenerates in one point
%
We_c = [Wnp, Wnp*gammaR_ref]; %add the torques to the regression matrix in order to perform the OLS
Sigma_c = (We_c'*We_c)/(n*N);
alpha_c = OLS(Sigma_c,par+1,-1);
%
%%
%%%%%%%%%%%%%%%%%%%% Perform the OLS to find the vertexes of the symplex in the parameter's space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rng(4+hh)
ntau = tn*randn(size(Wns,1),1);
tau =  Wnp*gammaR_ref;
taus = tau + ntau;
We = [Wns, taus]; %add the torques to the regression matrix in order to perform the OLS
% ONE SHOT ID
Sigma = (We'*We)/(n*N);
alpha = OLS(Sigma,par+1,-1)/scalingF;
[min(alpha(1:end-1,:),[],2),gammaR_ref,max(alpha(1:end-1,:),[],2)]

%
m = 2000;
Xi = We;%(1:m,:);
Si = (Xi'*Xi)/length(Xi);
a = OLS(Si,par+1,-1);
a_ls0 = ((Xi(:,1:end-1)'*Xi(:,1:end-1))^-1)*Xi(:,1:end-1)'*Xi(:,end);
iter = N-m+1;
l0 = min(a(1:end-1,:),[],2); l0(1:5) = [0;0;0;0;0];
u0 = max(a(1:end-1,:),[],2);
%
lb = zeros(length(l0),iter); lb(:,1) = l0;
ub = zeros(length(u0),iter); ub(:,1) = u0;
%
% RLS
P = 1e2*eye(par);
P2 = 1e2*eye(par);
a_rls= ones(par,iter); %a_rls(:,1) = a_ls0;
a_dprls= ones(par,iter);
for i = 2:iter
    Xi = We(i:i+m,:);
    Si = (Xi'*Xi)/m;
    a = OLS(Si,par+1,-1);
    lb(:,i) = max(lb(:,i-1),min(a(1:end-1,:),[],2));
    ub(:,i) = min(ub(:,i-1),max(a(1:end-1,:),[],2));
    idx = lb(:,i)>ub(:,i);
    if any(idx)
        %temp = lb(idx,i);
        %lb(idx,i) = ub(idx,i);
        %ub(idx,i) = temp;
        lb(:,i) = lb(:,i-1);
        ub(:,i) = ub(:,i-1);
    end
    % RLS Update
    gamma = Xi(end,1:end-1)';
    y = Xi(end,end);
    K = P*gamma*((1+(gamma')*P*(gamma))^-1);
    a_rls(:,i) = a_rls(:,i-1) + K*(y-(gamma')*a_rls(:,i-1));
    P=(eye(par)-K*gamma')*P;
    
    % DPRLS Update
    gamma2 = Xi(end,1:end-1)';
    y2 = Xi(end,end);
    K2 = P2*gamma2*((1+(gamma2')*P2*(gamma2))^-1);
    a_dprls(:,i) = a_dprls(:,i-1) + K2*(y2-(gamma2')*a_dprls(:,i-1));
    P2=(eye(par)-K2*gamma2')*P2;
    idx_l = a_dprls(:,i)<lb(:,i);
    idx_u = a_dprls(:,i)>ub(:,i);
    if any(idx_l)
        a_dprls(idx_l,i) = lb(idx_l,i);
    end
    if any(idx_u)
        a_dprls(idx_u,i) = ub(idx_u,i);
    end
end
%
figure(5)
for i = 1:9
    subplot(330+i)
        %{
        plot(lb(i,:),'LineWidth',2)
        hold on
        plot(ub(i,:),'LineWidth',2)
        plot(a_rls(i,:),':g','LineWidth',2)
        plot(a_dprls(i,:),':c','LineWidth',2)
        plot([0,iter],[gammaR_ref(i),gammaR_ref(i)],'--k','LineWidth',2)
        hold off
        %}
        %semilogx(a_rls(i,:),':g','LineWidth',2)
        semilogx(lb(i,:),'LineWidth',2)
        hold on
        
        semilogx(ub(i,:),'LineWidth',2)
        %semilogx(a_dprls(i,:),':c','LineWidth',2)
        semilogx([1,iter],[gammaR_ref(i),gammaR_ref(i)],'--k','LineWidth',1)
        hold off
        box on
        ylabel(['p_',num2str(i)])
end

options = optimoptions('fmincon',...
                        'Display','iter',...
                        'Algorithm', 'sqp',...
                        'MaxFunctionEvaluations', 1*1e3);
x0 = a_rls(:,end);
x = fmincon(@(x) selection_cost(x,Wns,taus),x0,[],[],[],[],lb(:,end),ub(:,end),[],options);
%
%
tau_rls = reshape(Wns*a_rls(:,end),n,length(q));
tau_x = reshape(Wns*x,n,length(q));
tauss = reshape(taus,n,length(q));

%
%
%
%
ex = tau_x'-tauss';
erls = tau_rls'-tauss';
norm_ex = sqrt(ex(:,1).^2 + ex(:,2).^2);
norm_ers = sqrt(erls(:,1).^2 + erls(:,2).^2);
figure(6)
%plot(tau_rls','r','LineWidth',2)
subplot(211)
plot(T,tauss','k','LineWidth',2)
hold on 
plot(T,tau_x','--b','LineWidth',1)
hold off
box on
subplot(212)
plot(T,norm_ex,'k','LineWidth',2)

disp('torque rec. error')
[selection_cost(x,Wns,taus)/norm(tauss),selection_cost(a_rls(:,end),Wns,taus)/norm(tauss)]
disp('par. error')
[norm(x-gammaR_ref)/norm(gammaR_ref),...
norm(a_rls(:,end)-gammaR_ref)/norm(gammaR_ref),...
norm(alpha(1:end-1,end)-gammaR_ref)/norm(gammaR_ref)]

function IDX = selection_cost(x,Wns,taus)
IDX = norm(Wns*x-taus) + norm(x);
end