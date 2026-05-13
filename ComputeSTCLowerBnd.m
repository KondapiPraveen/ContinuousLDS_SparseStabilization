% Compute the Lower Bound on Inter-Event Times for ST Control
%{
Inputs:
    System Matrices: A,B
    Feedback Matrix: K
    Ordering of rows of K: K_part
    Lyapunov Derivative Matrix: Q0 % A0^{\T}*P + P*A0 + Q0 = 0
    Sparsity Levels: S
    decay factor: beta \in (0,1)
%}
%{
Output: tminS : lower bound on Inter-Event Time for all S
%}

function tminS = ComputeSTCLowerBnd(A,B,K,P,K_part,Q0,S,beta)
    lenS = numel(S); meff = numel(K_part); n = size(A,1); tminS = zeros(lenS,1);
    Fs = zeros(2*n,2*n); Fs(1:n,1:n) = beta*(A+B*K); % Auxilary system

    for i=1:lenS
        s = S(i); % sparsity

        % Initialization
        L = ceil(meff/s); % # Subsystems
        Fl = zeros(2*n,2*n,L); % Dynamics of [x, x_e]
        Ql = zeros(n,n,L); %  Matrix in the derivative of V(t) 

        % Obtain optimal feedback matrices and cvx combinations
        [Kl,~] = OptimalFeedback(K,K_part,s); % Get Optimal Feedback and Cvx. Comb.
        
        tmin = inf;
        for l=1:L
            % Compute Fl
            Beff = B*Kl(:,:,l); Aeff = A+Beff;
            Fl(:,:,l) = [Aeff, Beff; -Aeff, -Beff];

            % Compute Ql = (A+B*K_l)'*P + P*(A+B*K_l)
            Ql(:,:,l) = Aeff'*P + P*Aeff;

            % Compute tminS,l
            Qs = (Ql(:,:,l) - beta*Q0);
            tminl = first_negative_eig_time(Fl(:,:,l),Fs,P,Qs);

            if tminl < tmin
                tmin = tminl;
            end
        end
        tminS(i) = tmin;
    end
end

function [t_star, lambda_star, info] = first_negative_eig_time(Fl,Fs,P,Qs,opts)
    % ================================================================
    % Computes
    %
    %   t_star = min { t > 0 : lambda_min(M(t)) < 0 }
    %
    % where
    %
    %   M(t) = N(t)' * (Ql - beta*Q) * N(t)
    %
    % INPUTS:
    % -------
    % Fl : Dynamics of [x,x_e] with lth feedback matrix 
    % Fs : Auxilary comparison system
    % P : Lyapunov matrix
    % Qs  : Ql-beta*Q - Matrix in the derivative of V(t)
    %
    % beta  : scalar
    %
    % opts  : structure with fields
    %         opts.t0        = initial time (default 0)
    %         opts.dt        = coarse search step
    %         opts.tmax      = maximum search time
    %         opts.tol       = bisection tolerance
    %
    % OUTPUTS:
    % --------
    % t_star       : first time when lambda_min < 0
    % lambda_star  : minimum eigenvalue at t_star
    % info         : diagnostic structure
    %
    % ================================================================
    
    % Default options
    if nargin < 5
        opts = struct;
    end
    
    if ~isfield(opts,'t0')
        opts.t0 = 3.4*1e-3;
    end
    
    if ~isfield(opts,'dt')
        opts.dt = 1e-3;
    end
    
    if ~isfield(opts,'tmax')
        opts.tmax = 1;
    end
    
    if ~isfield(opts,'tol')
        opts.tol = 1e-8;
    end
    
    if ~isfield(opts,'use_eigs')
        opts.use_eigs = false;
    end
    
    % Compute M(t)
    Mt = @(t) Mfun(Fl,Fs,P,t);

    % Precompute symmetric matrix
    Qs = (Qs + Qs')/2;

    % Define scalar function
    f = @(t) min_eig_value(t, Mt, Qs);
    
    % Step 1: Coarse forward scan to detect sign change
    
    t_left  = opts.t0;
    f_left  = f(t_left);
    
    found = false;
    
    while t_left < opts.tmax
    
        t_right = t_left + opts.dt;
        f_right = f(t_right);
    
        % Detect first sign crossing
        if (f_left >= 0) && (f_right < 0)
            found = true;
            break;
        end
    
        % Move interval forward
        t_left = t_right;
        f_left = f_right;
    end
    
    % If no crossing found
    if ~found
        warning('No negative eigenvalue detected up to tmax');
    
        t_star = opts.tmax;
        lambda_star = NaN;
    
        info.status = 'NOT_FOUND';
        return;
    end
    
    % Step 2: Bisection refinement
    %{
    while (t_right - t_left) > opts.tol
    
        t_mid = 0.5*(t_left + t_right);
        f_mid = f(t_mid);
    
        if f_mid < 0
            t_right = t_mid;
        else
            t_left = t_mid;
        end
    end
    %}
    
    % Final estimate
    
    t_star = t_left;
    if t_left == opts.t0
        warning('Lower Bound is t0. Try again with reduced the start time');
    end
    lambda_star = f(t_star);
    
    % Diagnostics
    info.status = 'SUCCESS';
    info.interval = [t_left, t_right];
    info.tolerance = opts.tol;
end

function lambda_min = min_eig_value(t, Mt, Qs)
    % Computes smallest eigenvalue of N(t)' * Qs * N(t)
    n = size(Qs,1);
    % Evaluate N(t)
    Mlt = (Mt(t)+Mt(t))/2;
    %{
    if eigs(Mlt,1,'smallestabs') > 2*1e-6
        N = zeros(n,1);
    else
        N = null(Mlt,2e-6);
        if isempty(N)
            N = zeros(n,1);
        end
    end
    %}

    if eigs(Mlt,1,'smallestabs') > 2*1e-6
        lambda_min = 0;
    else
        %
        N = null(Mlt,1e-3);
        if isempty(N)
            lambda_min = 0;
        else
            % Form symmetric matrix
            M = N' * Qs * N;
            M = (M + M')/2;
            
            %lambda_min = min(eig(M));
            lambda_min = eigs(M,1,'smallestreal');
        end
        
        %}
        %lambda_min = -1;
    end
end

function Mlt = Mfun(Fl,Fs,P,t)
    n = size(P,1);
    S = [eye(n), zeros(n)];
    SEFl = S*expm(Fl*t);
    SEFs = S*expm(Fs*t);
    Mlt = S*(SEFl'*P*SEFl - SEFs'*P*SEFs)*S';
end