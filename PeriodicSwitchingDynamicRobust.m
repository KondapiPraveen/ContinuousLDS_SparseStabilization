%{
Inputs:
    System Matrices : A,B,C
    Feedback Matrices : K
    Luenberger gain : F
    Initial State : x0
    Row ordering of feedback matrix : K_part
    Sparsity : S
    NoiseBound : b
    Max Time of Simulation : tmax
%}
% Output: rhoS: Switching Period for different sparsity s in S
function rhoS = PeriodicSwitchingDynamicRobust(A,B,C,K,F,K_part,x0,S,b,tmax)
    lenS = numel(S); meff = numel(K_part); n = size(A,1); rhoS = zeros(lenS,1);
    p=size(C,1); matI = eye(2*n); %matW = [eye(n); -eye(n)];
    fig1 = figure; fig2 = figure;
    Markers = {'s','d','o','+','x','*','^','v'}; % For plotting
    
    for i=1:lenS
        s = S(i); % sparsity
        
        % Compute rho: periodic switching time
        rho = 1; rho = ComputeRho(A,B,K,K_part,s);
        rhoS(i) = rho;

        % Initialization
        L = ceil(meff/s); % # Subsystems
        NumPeriods = ceil(tmax/rho); % Number of Periods
        Atl = zeros(2*n,2*n,L); % store transfer matrix of dynamics of [x;e]
        InvAtl = zeros(2*n,2*n,L); % Store inverse of InvAtl
        isSingular = 0; % Check for Atl is invertible
        EpAtl = zeros(2*n,2*n,L); % Store matrix exponentials of Atls
        Phi = eye(2*n,2*n); % State Transition Matrix for a Period
        x_curr = [x0;-x0]; xNorm2 = zeros(NumPeriods+1,1); xNorm2(1) = norm(x0);
        time = 0:NumPeriods; MarkerInd = 1:2:NumPeriods;

        % Obtain optimal feedback matrices and cvx combinations
        [Kl, OptCvx] = OptimalFeedback(K,K_part,s); % Get Optimal Feedback and Cvx. Comb.
        
        for l=1:L
            Beff = B*Kl(:,:,l);
            % Transfer matrix of the dynamics for [x;e]
            Atl(:,:,l) = [A+Beff, Beff; zeros(n,n), A + F*C];
            if rank(Atl(:,:,l))<2*n
               isSingular = 1;
            end
            EpAtl(:,:,l) = expm(Atl(:,:,l)*OptCvx(l)*rho);
            % State Transition for each period to compute the next state
            Phi = EpAtl(:,:,l)*Phi;
        end
        
        % Compute the state at the end of each period
        for k=1:NumPeriods
            x_curr = Phi*x_curr;
            xNorm2(k+1) = norm(x_curr(1:n));
        end
        
        xNorm2 = xNorm2.^2; % Get the sq. l2 norm of the state
        % Plotting
        figure(fig1)
        str = sprintf('$s = %d$',s);
        semilogy(time*rho,xNorm2,'DisplayName',str,'Marker',Markers(i),'MarkerIndices',MarkerInd,'MarkerSize',10,'LineWidth',2.5);
        hold on
        
        if isSingular == 1
            disp("The closed loop switching dynamics is singular. Use different approach for robustness plot");
            continue;
        end
        for l=1:L
            InvAtl(:,:,l) = Atl(:,:,l)\matI;
        end
        
        % Compute state at the end of each period under bounded noise
        xnNorm2 = zeros(NumPeriods+1,1); xnNorm2(1) = norm(x0);
        xn_inter = [x0;-x0]; % Intermediate state after each switching
        % Same noise vector in each period
        NoiseVec1 = b*rand(n,NumPeriods); % Process noise
        NoiseVec2 = b*rand(p,NumPeriods); % Measurement noise
        
        for k=1:NumPeriods
            NoiseNow = [NoiseVec1(:,k); -F*NoiseVec2(:,k)-NoiseVec1(:,k)];
            for l=1:L
                xn_inter = EpAtl(:,:,l)*xn_inter + InvAtl(:,:,l)*(EpAtl(:,:,l)-matI)*NoiseNow;
            end
            xnNorm2(k+1) = norm(xn_inter(1:n));
        end

        xnNorm2 = xnNorm2.^2; % Get the sq. l2 norm of the state
        % Plotting
        figure(fig2)
        str = sprintf('$s = %d$',s);
        semilogy(time*rho,xnNorm2,'DisplayName',str,'Marker',Markers(i),'MarkerIndices',MarkerInd,'MarkerSize',10,'LineWidth',2.5);
        hold on
    end
    figure(fig1);
    legend('NumColumns',2,'Interpreter','latex','FontSize',16); 
    xlabel('time (s)','Interpreter','latex','FontSize',16); ylabel('$\Vert \textbf{x}(t) \Vert^2$','Interpreter','latex','FontSize',16);
    title('Noiseless Dynamic Feedback')
    figure(fig2)
    legend('NumColumns',2,'Interpreter','latex','FontSize',16);
    xlabel('time (s)','Interpreter','latex','FontSize',16); ylabel('$\Vert \textbf{x}(t) \Vert^2$','Interpreter','latex','FontSize',16);
    title('Robust Dynamic Feedback')
end