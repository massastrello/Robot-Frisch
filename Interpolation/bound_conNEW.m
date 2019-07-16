function [c,ceq] = bound_conNEW(x,q,dq,ddq,B,ts,j_0)
   xi = abs(x(1))+ts;
   q = [q;x(2:3)];
   dq = [dq;x(4:5)];
   ddq = [ddq;x(6:7)];
   dddq = [j_0;x(8:9)];
   [qi, qdi, qddi] = eptapolytraj(xi,q',dq',ddq',dddq',ts);
    %{
    figure(4)
        subplot(311)
            plot(0:ts:xi,qi,'LineWidth',2)
            hold on
            plot([0,xi],[B(1),B(1)],'--k')
            plot([0,xi],[B(7),B(7)],'--k')
            hold off
        subplot(312)
            plot(0:ts:xi,qdi,'LineWidth',2)
            hold on
            plot([0,xi],[B(3),B(3)],'--k')
            plot([0,xi],[B(10),B(10)],'--k')
            hold off
        subplot(313)
            plot(0:ts:xi,qddi,'LineWidth',2)
            hold on
            plot([0,xi],[B(5),B(5)],'--k')
            plot([0,xi],[B(12),B(12)],'--k')
            hold off
    %}
    %
    c = [];
    %c = [max(max(abs(qi),[],2))-B(7),max(max(abs(qdi),[],2))-B(10),max(max(abs(qddi),[],2))-B(12)];
    c = [-min(qi,[],2)'+B(1:2),-min(qdi,[],2)'+B(3:4),-min(qddi,[],2)'+B(5:6),...
         max(qi,[],2)'-B(7:8), max(qdi,[],2)'-B(9:10),max(qddi,[],2)'-B(11:12)];
    ceq = [];
  
end