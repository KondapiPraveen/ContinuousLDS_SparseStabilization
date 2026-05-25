% Stability with 1-sparse input for several masses

% Choose the # masses
clear; clc; close all; rng(0);
tic;
Nl = 5:5:40; lenN = numel(Nl);
S = 1; tmax = 10; % Parameters of simulation
Markers = {'s','d','o','+','x','*','^','v'}; % For plotting
tminS = zeros(lenN,1); fig1 = figure;

for l=1:lenN
    N = Nl(l);
    c = zeros(N,1); n = 2*N;
    % Entries in the main, sub/super diagonal in symmetric toeplitz matrix
    %c(1) = -2; c(2) = 1;
    c(1) = -6; c(2) = 3;
    T = toeplitz(c);
    
    % System matrices generation: A = [0 I;T 0], B = [0;I], C = [0 I] 
    A = zeros(n,n); % Transfer matrix
    A(1:N,N+1:n) = eye(N);
    A(N+1:n,1:N) = T;
    
    B = zeros(n,N); % Input matrix
    % Effective Input Matrix
    %B_eff = eye(N);
    B_eff = diag(0.5*(1.5 + rand(N,1))); % Unif(0.75,1.25)
    B(N+1:n,:) = B_eff; m = N;
    
    p = N; C = zeros(p,n);
    C(:,N+1:n) = eye(N);
    
    % run PBH test to check if (A,B) is stabilizable
    % Ths test is recommended to run for new systems
    checkStable = 1;
    if checkStable
        isStable = isStabilizable(A,B);
        if isStable == 0
            error("The system is not stabilizable. Please choose a stabilizable system");
        end
    
        isDetect = isStabilizable(A.',C.');
        if isDetect == 0
            error("The system is not detectable. Please choose a detectable system");
        end
    end
    
    % Compute the stabilizing (non-sparse) controller using lqr MATLAB function
    % Find K s.t. A-BK is hurwitz using lqr function. Control input: u=-Kx
    Q = eye(n); R = eye(m); % state and input, cross regularization for state and input
    N = zeros(n,m); % cross regularization term PS : DO NOT CHANGE
    [K, P1] = lqr(A,B,Q,R,N); % Find the feedback matrix
    Ry = 5*eye(p); Ny = zeros(n,p);
    F = lqr(A.',C.',Q,Ry,Ny); % Find the Luenberger observer gain
    
    % Compute the lyapunov function V = x^{\T}Px s.t. (A-BK)^{\T}P + P(A-BK)
    % We compute this using lyap fuction in MATLAB
    K = -K; % In our work we consider A+BK is Hurwitz.
    A0 = A+B*K; % Closed loop dynamics
    Q0 = eye(n); % A0^{\T}*P + P*A0 + Q0 = 0
    P = lyap(A0',Q0); % lyap solves for X in A*X + X*A^{\T} + Q = 0, for A,Q inputs
    F = -F.'; % In our work we consider A+FC is Hurwitz
    A1 = A + F*C; % Error Dynamics
    PO = lyap(A1',Q0);
    
    % Compute Optimal Partition
    K_row_l2 = vecnorm(K,2,2).^2; % Compute l2 norm of each row of K
    [~,K_part] = sort(K_row_l2,'descend');
    meff = m - sum(K_row_l2 == 0); % number of all zero rows in K
    K_part_eff = K_part(1:meff);
    % State Feedback Continuous Control
    % Initialization
    
    %x0 = ones(n,1)/sqrt(n); % Initial state with l_2^2 2
    c4 = eigs(Q0,1,'smallestabs');
    assert(c4 > 0,'System is not stabilizable via switching') % Otherwise system is not stabilizable via switching
    gamma = c4/2; % gamma < c4
    assert(c4 > gamma,'Non Zeno behavior is not guaranteed'); % To avoid Zeno behavior
    % Comment to skip this simulation
    x0 = randn(n,1); x0 = x0/norm(x0); % Random vector on unit l1 ball
    [tminS(l),~,~,tl,xl,cntr] = StateFeedbackContinuousControl(A,B,K,P,x0,K_part_eff,S,gamma,tmax); close
    figure(fig1);
    MarkerInd = 1:10:cntr-1;
    semilogy(tl{1}(1:cntr-1),xl{1}(1:cntr-1),'LineWidth',2.5,'Marker',Markers(l),'MarkerSize',10,'MarkerIndices',MarkerInd); hold on;
end
legend([strcat('$N = $',num2str(Nl.'))],'interpreter','latex','FontSize',16,'NumColumns',2); grid on;
xlabel('time $(t)$','Interpreter','latex','FontSize',16); ylabel('$V(t)$','Interpreter','latex','FontSize',16);
toc;