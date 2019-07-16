%% Last Edit 26/04/2017
function [gammaF,IC1,IC2,T] = ModelReduction(q,qd,qdd,gamma,n,N)
global AI
global NI
global IC
global ControlCode
NI = [];
AI = [];
IC = [];
%% First SVD
disp('* Performing a first SVD of the regression matrix built upon an optimized trajectory for the entire model')
Wn = computeRegression(q(2:N-1,:),qd(2:N-1,:),qdd(2:N-1,:),n,N-2);
[~,S,V] = svd(Wn);
p = 0;
c = size(Wn,2);
for i = 1:c
    if norm(S(:,i))>1e-4
        p = p+1;
    end  
end
if p < c
    disp(['   The regression matrix is rank deficient ', num2str(p), '/', num2str(c),': model reduction needed'])
else
    disp(['   The regression matrix if full rank', num2str(p), '/', num2str(c),': no model reduction needed']);
    return
end
%
V1 = V(:, 1:p);
V1t = V1';
V2 = V(:, p+1:c);
V2t = V2';
% Find the Absolutely Identifiable Par.s 
for i = 1:c
    j = 1;
    Id = 0;
    while (j <= (c - p)) && (Id == 0)
        if abs(V2t(j,i)) > 1e-4
             Id = 1;
        else
            j = j + 1;
        end
    end
    if not(Id) 
        AI = [AI i];
    else
        IC = [IC,i];
    end
end
disp(['   ',num2str(size(AI,2)),' Absolutely Identifiable (AI) parameters found']);
% Find the Non Identifiable Par.s 
for i = 1:c
    j = 1;
    Id = 0;
    while (j <= p) && (Id == 0)
        x = abs(V1t(j,i));
        if x > 1e-4
             Id = 1;
        else
            j = j + 1;
        end
    end
    if not(Id) 
        NI = [NI i];
    end
end
disp(['   ',num2str(size(NI,2)),' NOT Identifiable (NI) parameters found']);
%% Rewrite the parameters vector
gammaA = gamma(AI);
gammaB = gamma;
gammaB(AI) = [];
%% Second SVD
disp('* Performing a second SVD of the regression matrix after deleting all the columns corresponding to the AI and NI parameters (Wr)')
Wr = ModcomputeRegression(q(2:N-1,:) , qd(2:N-1,:) , qdd(2:N-1,:));
[~,Sp,Vp] = svd(Wr);
c1 = size(Wr,2);
p1 = 0;
for i = 1:c1
    if norm(Sp(:,i))>1e-4
        p1 = p1+1;
    end  
end
% Check that p1 is equal to c - p
if not(p1==p-size(AI,2))
    disp(['   The rank of Wr is ' , num2str(p1), ', while the columns needed to complete the regression matrix are', num2str(p-size(AI,2))]);
else
    disp(['   The rank of Wr is ' , num2str(p1), ', equal to the columns needed to complete the regression matrix']);
end
V2p = Vp(:, p1+1:c1);

check = Wr*V2p;
if norm(check)> 1e-3
    error('   Norm of Wr*V2'' is not 0 (%f)', norm(check));
end
%% Replacing all the "numerical zeros" in V2p
disp('   Replacing all the "numerical zeros" of V2''...');
[V2p,count] = ReplaceZeros(V2p,1e-4);%%%%%%%OCCHIO!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%
disp('   Done.');
disp(['   ',num2str(count),'entries of V2'' substituted']);
%% Find Permutation Matrix
% Find all possible permutations: they are actually the n!/k!(n-k)! combinations of k
% el.s of an n size set.
disp(['   Retriving all the possible combinations of ',num2str(size(V2p,2)), ' rows of V2'' (over ',num2str(size(V2p,1)),')...']);
%mex -O VChooseK.c;
%Cn = VChooseK(uint8(1:size(V2p,1)),size(V2p,2));
%load Cn;
%disp('   Done.');
askopti = isequal(ControlCode,'ctrl');
if askopti
    opti = input('Do you like to optimize the permutation? [''YES''/''NO'']');
        disp('Optimizing for the best permutation ...');
        [bestC,best] = AllCombinations(size(V2p,1),size(V2p,2),V2p,opti);
        disp('   Done.');
        disp(['   The best cost function value found is ',num2str(best),' (at iteration ',num2str(i),')']);
else
    opti = 'YES';
    disp('Optimizing for the best permutation ...');
    [bestC,best] = AllCombinations(size(V2p,1),size(V2p,2),V2p,opti);
    disp('   Done.');
    disp(['   The best cost function value found is ',num2str(best),' (at iteration ',num2str(i),')']);
end
V22 = zeros(size(V2p,2),size(V2p,2));  
V11 = V2p;
gammaB1 = gammaB;
gammaB2 = sym(zeros(size(V2p,2),1));
IC1 = IC;
IC2 = zeros(1,size(V2p,2));
i = size(V2p,2);
while  i > 0
   V22(i,:) = V2p(bestC(i),:);
   V11(bestC(1,i),:) = [];
   gammaB2(i,:) = gammaB(bestC(i),:);
   gammaB1(bestC(1,i),:) = [];
   IC2(i) = IC(bestC(1,i));
   IC1(bestC(1,i)) = [];
   i = i - 1;
end
check = det(V22);
if check
    disp(['   CHECK: Optimal Permutation found! det(V22) = ',num2str(det(V22))]);
else
    error('   The optimal permutation found is actually bad: det(V22) = 0');
end
if not(isequal(IC,sort([IC1,IC2])))
    error('GodDamned!! error code: (-1)');
end

T = V11*(V22^-1);
T = ReplaceZeros(T,1e-2); %Replace numerical Zeros
gammaR = gammaB1 - sym(T)*gammaB2;
gammaF = [gammaA;gammaR];
digits 6;
vpa(gammaF);
disp('* Model Reduction Completed!!');
end
