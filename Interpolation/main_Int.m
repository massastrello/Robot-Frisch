
B = [-3,-3,-5,-5,-10,-10,3,3,5,5,10,10];
npts = size(q,1);
ts = 0.001;
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
options = ...%optimoptions('patternsearch','Display','iter');
optimoptions('fmincon',...
             'Display','iter',...
             'Algorithm', 'sqp',...
             'MaxFunctionEvaluations', 1*1e3);
                        % 'Algorithm', 'sqp',...
                        
i = 1;
th = 1e-2;
thi = th;
cnt = 0;
j_0 = [0,0];
tpts_opt = zeros(npts-1,1);
jpts_opt = zeros(npts-1,2);
while i < npts
    disp(['Iteration: ',num2str(i),'/',num2str(npts)])
    x0 = rand(1,3);
    [x,fval,flag] = fmincon(@(x) cost_int_time7(x),x0,[],[],[],[],[],[],...
        @(x) bound_con7(x,q(i:i+1,:),qd(i:i+1,:),qdd(i:i+1,:),B,ts,j_0),options);
    %[x,fval,flag] = patternsearch(@(x) cost_int_time7(x),x0,[],[],[],[],...
    %    [0.5,-inf,-inf],[inf,inf,inf],...
    %    @(x) bound_con7(x,q(i:i+1,:),qd(i:i+1,:),qdd(i:i+1,:),B,ts,j_0),options);
    %
    if (flag==-2)||(flag==0)
        disp('Unfeasible...')
        max_con = max(bound_con7(x,q(i:i+1,:),qd(i:i+1,:),qdd(i:i+1,:),B,ts,j_0));
        if max_con > thi
            disp(['max. constraint ABOVE threshold(',num2str(max_con),'/',num2str(thi),')'])
            q(i+1:end,:) = q([end,i+1:end-1],:);
            qd(i+1:end,:) = qd([end,i+1:end-1],:);
            qdd(i+1:end,:) = qdd([end,i+1:end-1],:);
            cnt = cnt + 1;
            disp('Attempting to re-assess sample')
            disp(['Attempt #:', num2str(cnt),'/',num2str((size(q,1)-i))])
            if (cnt-(size(q,1)-i))==0
                i = i - 1;
                thi = 2*th;
                cnt = 0;
                p = i + randperm(length(q(i+1:end,1)));
                q(i+1:end,:) = q(p,:);
                qd(i+1:end,:) = qd(p,:);
                qdd(i+1:end,:) = qdd(p,:);
            end
        else
            cnt = 0;
            cnt_flag = 0;
            ras_cnt = 0;
            %thi = th;
            i = i+1;
            disp(['max. constraint BELOW threshold(',num2str(max_con),'/',num2str(th),')'])
            tpts_opt(i) = x(1);
            jpts_opt(i,:) = x(2:3);
            j_0 = x(2:3);
        end
    else
        cnt = 0;
        cnt_flag = 0;
        ras_cnt = 0;
        %thi = th;
        i = i+1;
        tpts_opt(i) = x(1);
        jpts_opt(i+1,:) = x(2:3);
        j_0 = x(2:3);
    end
end


%Tpts = zeros(npts,1);
t0 = 0;
qint = [];
qdint = [];
qddint = [];
T = [];
Tpts = t0;
j = i;
for i= 1:j-1
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
%
save('traj_int.mat','qint','qdint','qddint','T','Tpts','q','qd','qdd')