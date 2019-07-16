function [y,count] = ReplaceZeros(x,tol)
count = 0;
for i = 1:size(x,1)
    for j = 1:size(x,2)
        if abs(x(i,j)) < tol
            x(i,j) = 0;
            count = count + 1;
        end
    end
end
y = x;
end
