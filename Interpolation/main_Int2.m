
B = [-3,-3,-5,-5,-10,-10,3,3,5,5,10,10];
npts = size(q,1);
ts = 0.01;
tpts = 1:size(q,1);

figure(3)
subplot(311)
    plot(tpts,q)
    hold on
    plot([0,npts],[B(1),B(1)],'--k')
    plot([0,npts],[B(7),B(7)],'--k')
    hold off
    box on
subplot(312)
    plot(tpts,qd)
    hold on
    plot([0,npts],[B(3),B(3)],'--k')
    plot([0,npts],[B(10),B(10)],'--k')
    hold off
    box on
subplot(313)
    plot(tpts,qdd)
    hold on
    plot([0,npts],[B(5),B(5)],'--k')
    plot([0,npts],[B(12),B(12)],'--k')
    hold off
    box on
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options1 = optimoptions('fmincon',...
                        'Algorithm', 'sqp',...
                        'MaxFunctionEvaluations', 1*1e4);
options2 =optimoptions('ga',...
                       'ConstraintTolerance',1e-6,...
                       'UseParallel', true,...
                       'Display', 'iter');optimoptions('patternsearch','Display','iter');

%optimoptions('ga',...
%                        'ConstraintTolerance',1e-6,...
%                        'UseParallel', true,...
%                        'Display', 'iter');
                        %'PlotFcn', @gaplotbestf);

%optimoptions('fmincon',...
           %             'Display','iter',...
           %             'MaxFunctionEvaluations', 1*1e4);
                        
                        % 'Algorithm', 'sqp',...
                        
i = 1;
th = 0.5;
cnt = 0;
j_0 = [0,0];
tpts_opt = zeros(npts-1,1);
jpts_opt = zeros(npts,2);
while i < npts
    disp(['Iteration: ',num2str(i),'/',num2str(npts)])
    x0 = rand(1,3);
    [x,fval,flag] = fmincon(@(x) cost_int_time7(x),x0,[],[],[],[],[],[],...
        @(x) bound_con7(x,q(i:i+1,:),qd(i:i+1,:),qdd(i:i+1,:),B,ts,j_0),options1); 
    %
    if (flag==-2)||(flag==0)
        x0 = rand(1,9);
        qi = q; qdi = qd; qddi = qdd;
        qi(i+1,:) = []; qdi(i+1,:) = []; qddi(i+1,:) = [];
        %[x,fval,flag] = fmincon(@(x) reassess_cost(x,qi,qdi,qddi),x0,[],[],[],[],...
        %    [0.05;B(1:6)';-100;-100],[100;B(7:12)';100;100],...
        %    @(x) bound_conNEW(x,q(i,:),qd(i,:),qdd(i,:),B,ts,j_0),options2);
        %rng default;
        [x,fval,flag] = ga(@(x) reassess_cost(x,qi,qdi,qddi),9,[],[],[],[],...
            [0.05;B(1:6)';-100;-100],[100;B(7:12)';100;100],...
            @(x) bound_conNEW(x,q(i,:),qd(i,:),qdd(i,:),B,ts,j_0),options2);
        %[x,fval,flag] = patternsearch(@(x) reassess_cost(x,qi,qdi,qddi),x0,[],[],[],[],...
        %    [0.05;B(1:6)';-100;-100],[100;B(7:12)';100;100],...
        %    @(x) bound_conNEW(x,q(i,:),qd(i,:),qdd(i,:),B,ts,j_0),options2);
        if (flag==-2)%||(flag==0)
            %{
            if i==1
                p = randperm(length(q));
                q = q(p,:);
                qd = qd(p,:);
                qdd = qdd(p,:);
            end
            %}
            i = i - 1;
            p = i + randperm(length(q(i+1:end,1)));
            q(i+1:end,:) = q(p,:);
            qd(i+1:end,:) = qd(p,:);
            qdd(i+1:end,:) = qdd(p,:);
            continue;
        else
            q(i+1,:) = x(2:3); qd(i+1,:) = x(4:5); qdd(i+1,:) = x(6:7);
            tpts_opt(i) = x(1);
            jpts_opt(i+1,:) = x(8:9);
            j_0 = x(8:9);
            i = i + 1; 
        end
    else
       tpts_opt(i) = x(1);
       jpts_opt(i+1,:) = x(2:3);
       j_0 = x(2:3);
       i = i + 1; 
    end
    %}
end


%Tpts = zeros(npts,1);
t0 = 0;
qint = [];
qdint = [];
qddint = [];
T = [];
Tpts = t0;
for i= 1:npts-1
    qddd = [jpts_opt(i,:);jpts_opt(i+1,:)]; 
    ti = abs(tpts_opt(i+1))+ts;
    [qint_i, qdint_i, qddint_i] = eptapolytraj(ti,q(i:i+1,:)',...
                                 qd(i:i+1,:)',qdd(i:i+1,:)',qddd',ts);
    
    qint = [qint,qint_i];
    qdint = [qdint,qdint_i];
    qddint = [qddint,qddint_i];
    T = [T,t0:ts:t0+ti];
    Tpts = [Tpts,t0+ti];
    t0 = t0 + ti;
    figure(4)
        subplot(311)
            plot(0:ts:ti,qint_i,'LineWidth',2)
            hold on
            plot([0,ti],[B(1),B(1)],'--k')
            plot([0,ti],[B(7),B(7)],'--k')
            hold off
        subplot(312)
            plot(0:ts:ti,qdint_i,'LineWidth',2)
            hold on
            plot([0,ti],[B(3),B(3)],'--k')
            plot([0,ti],[B(10),B(10)],'--k')
            hold off
        subplot(313)
            plot(0:ts:ti,qddint_i,'LineWidth',2)
            hold on
            plot([0,ti],[B(5),B(5)],'--k')
            plot([0,ti],[B(12),B(12)],'--k')
            hold off
        pause(2)
end
%
figure(4)
subplot(311)
    plot(Tpts,q,'*k')
    hold on
    plot(T,qint,'LineWidth',2)
    plot([0,Tpts(end)],[B(1),B(1)],'--k')
    plot([0,Tpts(end)],[B(7),B(7)],'--k')
    hold off
    box on
subplot(312)
    plot(Tpts,qd,'*k')
    hold on
    plot(T,qdint,'LineWidth',2)
    plot([0,Tpts(end)],[B(3),B(3)],'--k')
    plot([0,Tpts(end)],[B(10),B(10)],'--k')
    hold off
    box on
subplot(313)
    plot(Tpts,qdd,'*k')
    hold on
    plot(T,qddint,'LineWidth',2)
    plot([0,Tpts(end)],[B(5),B(5)],'--k')
    plot([0,Tpts(end)],[B(12),B(12)],'--k')
    hold off
    box on
%
save('traj_int.mat','qint','qdint','qddint','T','Tpts','q','qd','qdd')