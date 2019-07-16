function IDX = cost_int_time(x,q,dq,ddq,B,ts)
    xi = abs(x)+ts;
    [qi, qdi, qddi, ~] = quinticpolytraj(q',[0,xi], 0:ts:xi,...
                            'VelocityBoundaryCondition', dq',...
                            'AccelerationBoundaryCondition', ddq');
    c = [max(max(abs(qi),[],2))-B(7),max(max(abs(qdi),[],2))-B(10),max(max(abs(qddi),[],2))-B(12)];
    IDX = abs(x) + sum(c); 
end