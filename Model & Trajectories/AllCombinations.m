function [best,fval]= AllCombinations(N,K,V2,opti)
m = 0;
h = K;
iteration = 1;         
A = 1:K;          
fval = 1e20;
tot = factorial(N)/(factorial(K)*factorial(N-K));
while(iteration>0)
    [A,m,h,iteration] = Combination(N, K, A, m, h, iteration);
    f = PermutationCost1(V2,A);
    if f<fval
        fval = f;
        best = A;
        if isequal(opti,'NO')
            break;
        end
    end
    if (iteration*100/tot == 0) || (iteration*100/tot == 25) || (iteration*100/tot == 50) || (iteration*100/tot == 75) || (iteration*100/tot == 100)
        disp(['   ',num2str(iteration*100/tot),'%']);
    end
end
end