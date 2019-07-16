function P = generateparticles(alpha,Np)
lb = min(alpha,[],2);
ub = max(alpha,[],2);
%
lb = lb(1:end-1);
ub = ub(1:end-1);
%
X = zeros(length(lb),Np);
for i = 1:length(lb)
    X(i,:) = linspace(lb(i),ub(i),Np);
end
%
counter = 1;
P = zeros(Np^(length(lb)),length(lb));
for a = 1:Np
    for b = 1:Np
        for c = 1:Np
            for d = 1:Np
                for e = 1:Np
                    for f = 1:Np
                        for g = 1:Np
                            for h = 1:Np
                                for i = 1:Np
                                    P(counter,:) = [X(1,a),X(2,b),X(3,c),...
                                                    X(4,d),X(5,e),X(6,f),...
                                                    X(7,g),X(8,h),X(9,i)];
                                    counter = counter + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end