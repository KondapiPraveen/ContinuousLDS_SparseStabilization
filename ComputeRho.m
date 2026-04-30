% Compute the periodic switching time (rho) of given system
%{
Inputs:
    System Matrices : A, B
    Feedback Matrix : K
    Row_Priority : K_part (rows with highest to lowest l2 norms)
    s : sparsity
%}
% Output: rho (periodic switching time)
function rho = ComputeRho(A,B,K,K_part,s)
    [Kl, OptCvx] = OptimalFeedback(K,K_part,s);
    meff = numel(K_part); L = ceil(meff/s); tol = 1e-6;
    
    n = size(A,1); Al = zeros(n,n,L);
    SumNorm = 0; % Sum of induced l2 norm of switched systems
    for l=1:L
        Al(:,:,l) = A + B*Kl(:,:,l);
        SumNorm = SumNorm + OptCvx(l)*norm(Al(:,:,l));
    end

    rho_prev = 1; rho_curr = log(2)/SumNorm; CheckStability = 1;
    % Check if the periodic switching is stable with rho_curr
    while CheckStability
        Phi_Rho = eye(n);
        % Compute the product of matrix exponentials
        for l=1:L
            Phi_Rho = expm(Al(:,:,l)*OptCvx(l)*rho_curr)*Phi_Rho;
        end
        isStable = norm(Phi_Rho) < 1;
        isFiner = (abs(rho_prev - rho_curr) < tol);
        CheckStability = (~isStable) || (~isFiner);
        if CheckStability
            rho_next = ((rho_curr+rho_prev)/2)*isStable + rho_curr/2*(~isStable);
            rho_prev = rho_prev*isStable + rho_curr*(~isStable);
            rho_curr = rho_next;
        end
    end
    rho = rho_curr;
end