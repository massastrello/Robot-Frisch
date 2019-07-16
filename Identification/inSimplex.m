function P_1 = inSimplex(P_0,alpha)
n = size(alpha,1)-1;
in = inhull(P_0,alpha(1:n,:)');
P_1 = in.*P_0;
P_1( ~any(P_1,2), : ) = [];
end