function [qi,qdi,qddi] = eptapolytraj(dt,q,qd,qdd,qddd,ts)
    n = size(q,1);
    T = 0:ts:dt;
    qi = zeros(n,size(T,2));
    qdi = zeros(n,size(T,2));
    qddi = zeros(n,size(T,2));
    for i = 1:n
        h = q(i,2)-q(i,1);
        a0 = q(i,1);
        a1 = qd(i,1);
        a2 = qdd(i,1)/2;
        a3 = qddd(i,1)/6;
        a4 = (210*h -...
              dt*((30*qdd(i,1)-15*qdd(i,2))*dt + ...
              (4*qddd(i,1)+qddd(i,2))*(dt^2)+...
              120*qd(i,1)+90*qd(i,2) ))/(6*(dt^4));
        a5 = (-168*h +...
              dt*((20*qdd(i,1)-14*qdd(i,2))*dt + ...
              (2*qddd(i,1)+qddd(i,2))*(dt^2)+...
              90*qd(i,1)+78*qd(i,2) ))/(2*(dt^5));
        a6 = (420*h -...
              dt*((45*qdd(i,1)-39*qdd(i,2))*dt + ...
              (4*qddd(i,1)+3*qddd(i,2))*(dt^2)+...
              216*qd(i,1)+204*qd(i,2) ))/(6*(dt^6));
        a7 = (-120*h +...
              dt*((12*qdd(i,1)-12*qdd(i,2))*dt + ...
              (qddd(i,1)+qddd(i,2))*(dt^2)+...
              60*qd(i,1)+60*qd(i,2) ))/(6*(dt^7));
       qi(i,:) = a0+a1.*T+a2.*(T.^2)+a3.*(T.^3)...
           +a4.*(T.^4)+a5.*(T.^5)+a6.*(T.^6)+a7.*(T.^7);
       qdi(i,:) = a1 + 2*a2.*T + 3*a3.*(T.^2) + 4*a4.*(T.^3)...
           +5*a5.*(T.^4) + 6*a6.*(T.^5) + 7*a7.*(T.^6);
       qddi(i,:) = 2*a2 + 6*a3.*T + 12*a4.*(T.^2) + 20*a5.*(T.^3)...
           +30*a6.*(T.^4) + 42*a7.*(T.^5);
    end
end