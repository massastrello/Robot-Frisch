% Computes the max. noise variance of each parameter from the matrix Sigma
function s = maxvar(Sigma)
    i = 1;
    s = zeros(size(Sigma,1),1);
    while (i<size(Sigma,1)+1)
    % Sigma_i is a matrix obtained deleting the i-th row and column of Sigma
        Sigma_i = Sigma;
        Sigma_i(i,:) = [];
        Sigma_i(:,i) = [];
        s(i,1) = det(Sigma)/det(Sigma_i); 
        i = i + 1;
    end
end