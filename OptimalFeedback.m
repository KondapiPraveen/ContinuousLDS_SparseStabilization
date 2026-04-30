% Compute Optimal Sparse Feedback Matrices
% Find out optimal cvx combination and optimal partition
function [Kl, OptCvx] = OptimalFeedback(K,K_part,s)
    meff = numel(K_part); L = ceil(meff/s);
    m = size(K,1); n = size(2);
    Kl = zeros(m,n,L); OptCvx = zeros(L);
    CONST = 0; % Find the optimal combination constant
    
    for l=1:L
        Sl = K_part((l-1)*s+1:min(l*s,meff));
        CONST = CONST + sqrt(sum(K_row_l2(Sl)));
    end
    
    for l = 1:L
        Sl = K_part((l-1)*s+1:min(l*s,meff));
        K_eff = zeros(m,n);
        K_eff(Sl,:) = K(Sl,:);
        %tol = max(abs(K_eff),[],"all");
        %K_eff(abs(K_eff) < tol*1e-3) = 0;
        OptCvx(l) = CONST/norm(K_eff,'fro');
        K_eff = OptCvx(l)*K_eff;
        Kl(:,:,l) = K_eff;
    end
end