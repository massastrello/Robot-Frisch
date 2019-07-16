%[A,M,H,ITERATION] = COMBINATION(N, K, A, M, H, ITERATION) returns 
%the next K-subset of an N-set. It implements the "Revolving Door" (NEXKSB)
%algorithm found in "Combinatorial Algorithms", 2nd Ed., by Albert Nijenhuis
%and Herbert S. Wilf. It is a recursive alternative to the NCHOOSEK command
%that is extremely fast (O(1)), at the expense of slightly more user
%interaction. (You do get something that NCHOOSEK does not give you for this
%extra effort: see input "m", below.) Apart from its efficiency, it is
%awesome because each new K-subset differs from the last by a single element
%and the last combination generated differs from the first combination by a
%single element.cdx



function [A,m,h,iteration,out] = Combination(N, K, A, m, h, iteration)
if(iteration~=1)
    h = (m>=(N-h))*h;
    h = h + 1;    
    m = A(K+1-h);
else
    m = 0;
    h = K;
end
A(K+(1:h)-h) = m + (1:h);
iteration = (A(1)~=(N-K+1))*(iteration + 1)-(A(1)==(N-K+1))*iteration;
end