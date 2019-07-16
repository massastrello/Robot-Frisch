% It takes as input the covariance matrix Sigma and gives back a matrix
% whose columns are the solutions of the Ordinary Least Sqares

function alpha = OLS(Sigma, norm, norm_val) %"norm" is the index of the variable that we want to normalize to norm_val (if norm = 0 & norm_val = 0--> no normalization)
i = 1;
alpha = zeros(size(Sigma,1),size(Sigma,1));
s = maxvar(Sigma);
while i < (size(Sigma,1)+ 1)
    a = zeros(size(Sigma,1),1);
    SigmaT = zeros(size(Sigma,1),size(Sigma,1));
    SigmaT(i,i) = s(i);
    SigmaH = Sigma - SigmaT;
    %corank = size(SigmaH,1) - rank(SigmaH);
        %if corank == 1 
            a = null(SigmaH);
            % Now normalize each set of parameters wrt one of the variables, if requested
            if norm > 0
                a = a*(norm_val/a(norm));
            end
        %else
        %    warning('corank(SigmaH) > 1 !!!');
        %end
    alpha(:,i) = a;
    i = i + 1;
end
end