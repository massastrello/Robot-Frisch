npts = size(q,1);
ts = 0.02;
scalef = 1;
Tpts = scalef.*(0:size(q,1)-1);
T = scalef.*(0:ts:size(q,1)-1);
[qint,qdint,qddint,pp] = quinticpolytraj(q',Tpts,T,...
                        'VelocityBoundaryCondition',qd',...
                        'AccelerationBoundaryCondition',qdd');

figure(4)
subplot(311)
    
    plot(T,qint,'LineWidth',2)
    hold on
    plot(Tpts,q,'*k')
    %xlim([0,30])
    %plot([0,Tpts(end)],[B(1),B(1)],'--k')
    %plot([0,Tpts(end)],[B(7),B(7)],'--k')
    hold off
    box on
subplot(312)
    plot(T,qdint,'LineWidth',2)
    hold on
    plot(Tpts,qd,'*k')
    %plot([0,Tpts(end)],[B(3),B(3)],'--k')
    %plot([0,Tpts(end)],[B(10),B(10)],'--k')
    %xlim([0,30])
    hold off
    box on
subplot(313)
    plot(T,qddint,'LineWidth',2)
    hold on
    plot(Tpts,qdd,'*k')
    xlim([0,30])
    %plot([0,Tpts(end)],[B(5),B(5)],'--k')
    %plot([0,Tpts(end)],[B(12),B(12)],'--k')
    hold off
    box on

%{
options = optimoptions('fmincon',...
             'Display','iter',...
             'MaxFunctionEvaluations', 1*1e3);
         %%'Algorithm', 'sqp',...

i = 1;
th = 1e-3;
thi = th;
cnt = 0;
j_0 = [0,0];
tpts_opt = zeros(npts-1,1);
%jpts_opt = zeros(npts-1,2);
jpts_opt = zeros(npts-1,4);
while i < npts
    disp(['Iteration: ',num2str(i),'/',num2str(npts)])
    x0 = rand(1,5);%[1,0,0];
    [x,fval,flag] = fmincon(@(x) cost_int_time7(x),x0,[],[],[],[],[],[],...
        @(x) bound_jerk(x,q(i:i+1,:),qd(i:i+1,:),qdd(i:i+1,:),B,ts,j_0),options);
    if (flag==-2)||(flag==0)
        disp('Unfeasible...')
        max_con = max(bound_con7(x,q(i:i+1,:),qd(i:i+1,:),qdd(i:i+1,:),B,ts,j_0));
        %in = input('Is this acceptable? (y/n)');
        if max_con > thi
            disp(['max. constraint ABOVE threshold(',num2str(max_con),'/',num2str(thi),')'])
            q(i+1:end,:) = q([end,i+1:end-1],:);
            qd(i+1:end,:) = qd([end,i+1:end-1],:);
            qdd(i+1:end,:) = qdd([end,i+1:end-1],:);
            cnt = cnt + 1;
            disp(['Attempting to re-assess sample'])
            disp(['Attempt #:', num2str(cnt),'/',num2str((size(q,1)-i))])
            if (cnt-(size(q,1)-i))==0
                cnt = 0;
                thi = 2*thi;
            end                
        else
            cnt = 0;
            cnt_flag = 0;
            ras_cnt = 0;
            i = i+1;
            disp(['max. constraint BELOW threshold(',num2str(max_con),'/',num2str(th),')'])
            tpts_opt(i) = x(1);
            %jpts_opt(i,:) = x(2:3);
            jpts_opt(i,:) = x(2:5);
            j_0 = x(2:3);
            thi = th;
        end
    else
        cnt = 0;
        cnt_flag = 0;
        ras_cnt = 0;
        i = i+1;
        tpts_opt(i) = x(1);
        %jpts_opt(i+1,:) = x(2:3);
        jpts_opt(i,:) = x(2:5);
        j_0 = x(2:3);
        thi = th;
    end
end

t0 = 0;
qint = [];
qdint = [];
qddint = [];
T = [];
Tpts = t0;
j = i;
for i= 1:j-1
    %qddd = [jpts_opt(i,:);jpts_opt(i+1,:)];
    qddd = [jpts_opt(i,2:3);jpts_opt(i,4:5)];
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
    plot(Tpts,q(1:j,:),'*k')
    hold on
    plot(T,qint,'LineWidth',2)
    plot([0,Tpts(end)],[B(1),B(1)],'--k')
    plot([0,Tpts(end)],[B(7),B(7)],'--k')
    hold off
    box on
subplot(312)
    plot(Tpts,qd(1:j,:),'*k')
    hold on
    plot(T,qdint,'LineWidth',2)
    plot([0,Tpts(end)],[B(3),B(3)],'--k')
    plot([0,Tpts(end)],[B(10),B(10)],'--k')
    hold off
    box on
subplot(313)
    plot(Tpts,qdd(1:j,:),'*k')
    hold on
    plot(T,qddint,'LineWidth',2)
    plot([0,Tpts(end)],[B(5),B(5)],'--k')
    plot([0,Tpts(end)],[B(12),B(12)],'--k')
    hold off
    box on
%}
%save('traj_int.mat','qint','qdint','qddint','T','Tpts','q','qd','qdd')
