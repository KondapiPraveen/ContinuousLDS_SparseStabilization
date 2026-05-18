% Compute the upper bound on the norm of state transition matrix for
% periodic switching

%{
Inputs: 
    Closed loop system matrices: Al = A + B*Kl
    Optimal Cvx Combination for each subsystem: OptCvx
    Switching Period: rho
%}

function [alpha,beta,k1] = ComputePhiNorm(Al,OptCvx,rho)
    n = size(Al,1); % size of the state
    L = size(Al,3); % # subsystems
    Del = 1e-3; % Decrease Del if L is large and rho is small
    if Del > min(OptCvx*rho)
        error("Choose smaller Del or increase $\rho$")
    end
    EpAl = zeros(n,n,L); EpAlDel = zeros(n,n,L);
    SwitchTimes = cumsum(OptCvx)*rho;
    
    % Compute Phi
    Phi = eye(n);
    for l=1:L
        EpAl(:,:,l) = expm(Al(:,:,l)*OptCvx(l)*rho);
        EpAlDel(:,:,l) = expm(Al(:,:,l)*Del);
        % State Transition for each period to compute the next state
        Phi = EpAl(:,:,l)*Phi;
    end

    Abar = logm(Phi);
    Q = eye(n);
    P = lyap(Abar',Q);
    Pmax = norm(P); Pmin = eigs(P,1,'smallestreal');
    % ||Phi(t)||^2 \leq beta*exp(-alpha*t)
    alpha = 1/Pmax; beta = Pmax/Pmin;

    % Compute k1 = \max_{0<=t<=rho} ||Phi(0,t)||
    t = Del; Phit = eye(n); k1 = 1; l_curr = 1;
    while t <= rho
        % search the completed subsytems
        l_prev = l_curr;
        BinarySwitch = (t <= SwitchTimes);
        l_curr = find(BinarySwitch,1);

        % Compute Phi(t)
        if l_prev ~= l_curr
            Phit = eye(n);
            for l=1:l_prev
                Phit = EpAl(:,:,l)*Phit;
            end
            Phit = expm(Al(:,:,l)*(t-SwitchTimes(l_prev)))*Phi;
        else
            Phit = EpAlDel(:,:,l_curr)*Phit;
        end
        t = t + Del; % Update Time
        nPhit = norm(Phit);
        if k1 < nPhit
            k1 = nPhit;
        end
    end
end