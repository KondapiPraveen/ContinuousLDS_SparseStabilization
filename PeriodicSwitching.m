% Periodic Switching under state feedback and compare with trajectory based
% switching
%{
Inputs:
    System Matrices : A,B
    Feedback Matrices : K
    Initial State : x0
    Row ordering of feedback matrix : K_part
    Sparsity : S
%}
% Output: rhoS: Switching Period for different sparsity s in S
function rhoS = PeriodicSwitching(A,B,K,x0,K_part,S)
    lenS = numel(S); meff = numel(K_part); n = size(A,1); figure();
    Markers = {'s','d','o','+','x','*','^','v'}; % For plotting
    legend('NumColumns',2,'Interpreter','latex');
    rhoS = zeros(lenS,1);
    for i=1:lenS
        s = S(i); % sparsity
        
        % Compute rho: periodic switching time
        % Improve time by precomputing rho
        rho = ComputeRho(A,B,K,K_part,s); rhoS(i) = rho;

        % Initialization
        L = ceil(meff/s); tmax = 10; % # Subsystems, Max Time of Simulation
        NumPeriods = ceil(tmax/rho); % Number of Periods
        Atl = zeros(n,n,L); % store transfer matrix of dynamics of [x;e]
        Phi = eye(n,n); % State Transition Matrix for a Period
        x_curr = x0; xNorm2 = zeros(NumPeriods+1,1); xNorm2(1) = norm(x0);

        % Obtain optimal feedback matrices and cvx combinations
        [Kl, OptCvx] = OptimalFeedback(K,K_part,s); % Get Optimal Feedback and Cvx. Comb.
        
        for l=1:L
            Beff = B*Kl(:,:,l);
            % Transfer matrix of the dynamics for [x;e]
            Atl(:,:,l) = A+Beff;
            % State Transition for each period to compute the next state
            Phi = expm(Atl(:,:,l)*OptCvx(l)*rho)*Phi;
        end
        
        % Compute the state at the end of each period
        for k=1:NumPeriods
            x_curr = Phi*x_curr;
            xNorm2(k) = norm(x_curr);
        end
        
        xNorm2 = xNorm2.^2; % Get the sq. l2 norm of the state
        % Plotting
        str = sprintf('$s = %d$',s);
        plot((NumPeriods-1)*rho,xNorm2,'DisplayName',str,'Marker',Markers(i),'MarkerSize',10,'LineWidth',2.5);
        hold on;
    end
end