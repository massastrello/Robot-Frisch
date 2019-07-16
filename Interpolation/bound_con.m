function [c,ceq] = bound_con(x,q,dq,ddq,B,ts)
    x = abs(x)+ts;
    [qi, qdi, qddi, ~] = quinticpolytraj(q',[0,x], 0:ts:x,...
                            'VelocityBoundaryCondition', dq',...
                            'AccelerationBoundaryCondition', ddq');
    %{
    figure(4)
        subplot(311)
            plot(0:ts:x,qi)
            hold on
            plot([0,x],[B(1),B(1)],'--k')
            plot([0,x],[B(7),B(7)],'--k')
            hold off
        subplot(312)
            plot(0:ts:x,qdi)
            hold on
            plot([0,x],[B(3),B(3)],'--k')
            plot([0,x],[B(10),B(10)],'--k')
            hold off
        subplot(313)
            plot(0:ts:x,qddi)
            hold on
            plot([0,x],[B(5),B(5)],'--k')
            plot([0,x],[B(12),B(12)],'--k')
            hold off
                          %}
    %
    c = [];
    c = [max(max(abs(qi),[],2))-B(7),max(max(abs(qdi),[],2))-B(10),max(max(abs(qddi),[],2))-B(12)];
    %c = [-min(qi,[],2)'+B(1:2),-min(qdi,[],2)'+B(3:4),-min(qddi,[],2)'+B(5:6),...
        % max(qi,[],2)'-B(7:8), max(qdi,[],2)'-B(9:10),max(qddi,[],2)'-B(11:12)];
    ceq = [];
  
end