function IDX = PermutationCost1(V2,Xi)
V22 = V2(Xi,:);
if abs(det(V22)) < 1e-6
    IDX = 1e20;
else
    IDX = max(abs(eig(V22)))/min(abs(eig(V22))) + 1/abs(det(V22));  
end
end