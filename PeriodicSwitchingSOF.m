%{
dx/st = Ax + Bu + v;  y = Cx + w
%}
%{
Inputs:
    System Matrices: A,B,C
    SOF Matrix: F
    row ordering of feedback matrix F*C : K_part
    Initial State: x0
    Sparsity : S
    NoiseBound : b
%}
%{
Output: rhoS: Switching Period for different sparsity s in S
%}

function rhoS = PeriodicSwitchingSOF(A,B,C,F,K_part,x0,S,b,tmax)
    K = F*C; % Feedback Matrix
    lenS = numel(S); meff = numel(K_part); n = size(A,1); rhoS = zeros(lenS,1);
    p=size(C,1); matI = eye(n);
    fig1 = figure;
    fig2 = figure;
    Markers = {'s','d','o','+','x','*','^','v'}; % For plotting

    for i=1:lenS
        s = S(i); % sparsity
        
        % Compute rho: periodic switching time
        rho = 1; rho = ComputeRho(A,B,K,K_part,s);
        rhoS(i) = rho;

        % Initialization
        L = ceil(meff/s); % # Subsystems
        NumPeriods = ceil(tmax/rho); % Number of Periods
        Atl = zeros(n,n,L); % store transfer matrix of dynamics of [x;e]
        InvAtl = zeros(n,n,L); % Store inverse of InvAtl
        Fl = zeros([size(F),L]); % Store Optimal SOF matrices
        BFl = zeros(n,p,L); % Store B*Fl
        isSingular = 0; % Check for Atl is invertible
        EpAtl = zeros(n,n,L); % Store matrix exponentials of Atls
        Phi = eye(n,n); % State Transition Matrix for a Period
        xNorm2 = zeros(NumPeriods+1,1); xNorm2(1) = norm(x0);
        time = 0:NumPeriods; MarkerInd = 1:2:NumPeriods;

        % Obtain optimal feedback matrices and cvx combinations
        [Kl, OptCvx] = OptimalFeedback(K,K_part,s); % Get Optimal Feedback and Cvx. Comb.

        for l=1:L
            Beff = B*Kl(:,:,l);
            % Transfer matrix of the dynamics for [x;e]
            Atl(:,:,l) = A+Beff;
            if rank(Atl(:,:,l))<n
               isSingular = 1;
            end
            EpAtl(:,:,l) = expm(Atl(:,:,l)*OptCvx(l)*rho);
            % State Transition for each period to compute the next state
            Phi = EpAtl(:,:,l)*Phi;
        end

        % Compute the state at the end of each period
        x_curr = x0;
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
            Sl = K_part((l-1)*s+1:min(l*s,meff));
            Fl(Sl,:,l) = OptCvx(l)*F(Sl,:);
            BFl(:,:,l) = B*Fl(:,:,l);
        end

        % Compute state at the end of each period under bounded noise
        xnNorm2 = zeros(NumPeriods+1,1); xnNorm2(1) = norm(x0);
        xn_inter = x0; % Intermediate state after each switching
        NoiseVec1 = b*rand(n,NumPeriods); % Process noise
        NoiseVec2 = b*rand(p,NumPeriods); % Measurement noise

        for k=1:NumPeriods
            for l=1:L
                NoiseNow = BFl(:,:,l)*NoiseVec2(:,k) + NoiseVec1(:,k);
                xn_inter = EpAtl(:,:,l)*xn_inter + InvAtl(:,:,l)*(EpAtl(:,:,l)-matI)*NoiseNow;
            end
            xnNorm2(k+1) = norm(xn_inter(1:n));
        end

        xnNorm2 = xnNorm2.^2; % Get the sq. l2 norm of the state
        % Plotting
        figure(fig2)
        str = sprintf('$s = %d$',s);
        plot(time*rho,xnNorm2,'DisplayName',str,'Marker',Markers(i),'MarkerIndices',MarkerInd,'MarkerSize',10,'LineWidth',2.5);
        hold on
    end
    figure(fig1)
    legend('NumColumns',2,'Interpreter','latex','FontSize',16); 
    xlabel('time (s)','Interpreter','latex','FontSize',16); ylabel('$\Vert \textbf{x}(t) \Vert^2$','Interpreter','latex','FontSize',16);
    title('Noiseless SOF')
    figure(fig2)
    legend('NumColumns',2,'Interpreter','latex','FontSize',16);
    xlabel('time (s)','Interpreter','latex','FontSize',16); ylabel('$\Vert \textbf{x}(t) \Vert^2$','Interpreter','latex','FontSize',16);
    title('Robust SOF')
end