% Plot Wrapper script
% CT-LDS : N mass-spring system 
% Find out more at http://www.ece.umn.edu/~mihailo/software/lqrsp/mass_spring.html
%
% In this script F is for Luenberger Observer Matrix
% In simulation under SOF, F is used for SOF matrix in the function defn 

close all; clear; clc;
% N denotes number of masses in the mass-spring system
N = 10; c = zeros(N,1); n = 2*N;
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
%% State Feedback Continuous Control
% Initialization
tic;
tmax = 10; % 10sec of simulation time
S = union(1,2:m/5:7);
x0 = ones(n,1)/sqrt(n); % Initial state with l_2^2 2
c4 = eigs(Q0,1,'smallestabs');
assert(c4 > 0,'System is not stabilizable via switching') % Otherwise system is not stabilizable via switching
gamma = c4/2; % gamma < c4
assert(c4 > gamma,'Non Zeno behavior is not guaranteed'); % To avoid Zeno behavior
% Comment to skip this simulation
tminS = StateFeedbackContinuousControl(A,B,K,P,x0,K_part_eff,S,gamma,tmax);
toc;
%% Feedback Gain under different partitions and convex combinations
tic;
% Choose Random B_eff
S = 1:1:m; % Sparsity levels
% Comment to skip this simulation
CtrlEnergyBound(K,K_part_eff,S);
toc;
%% Periodic Switching for Dynamic w/0 and w/ Bounded Noise
tic;
S = union(1,2:m/5:4*m/5);
b = 1/n; % Bound on each noise entry
x0 = ones(n,1)/sqrt(n); % Initial state with l_2^2 2
tmax = 20; % 500 sec simulation
% Comment to skip the simulation
rhoSDynamic = PeriodicSwitchingDynamicRobust(A,B,C,K,F,K_part,x0,S,b,tmax);
toc;
%% Periodic Switching for SOF w/0, w/ Bounded Noise
tic;
S = union(1,2:m/5:4*m/5); % Sparsity levels
x0 = ones(n,1)/sqrt(n); % Initial state with l_2^2 2
b = 1/n; % Bound on each noise entry
% Compute the SOF L
P1 = P1/2; maxItr = 100; thrsh = 0.1;
L = SOF(A,B,C,Q,R,N,P1,maxItr,thrsh);
L = -L; % In our case A+B*L*C is stable
Keff = L*C;
% Compute L_part
L_row_l2 = vecnorm(Keff,2,2).^2; % Compute l2 norm of each row of K
[~,L_part] = sort(L_row_l2,'descend');
% Comment to skip this simulation
tmax = 30; % 30sec of simulation time
rhoSSOF = PeriodicSwitchingSOF(A,B,C,L,L_part,x0,S,b,tmax);
toc;
%% Self Triggered Stabilization with State Feedback
tic;
S = union(1,2:m/5:4*m/5); % Sparsity levels
%S = m;
beta = 0.5; % Decay factor such that inter-event time > 0
tminSST = ComputeSTCLowerBnd(A,B,K,P,K_part,Q0,S,beta);
tic;
%% Simulation
tic;
x0 = ones(n,1)/sqrt(n); % Initial state with l_2^2 2
tmax = 10; % 10sec of simulation time
rmax = 10;
STSparseStabilization(A,B,K,K_part,P,x0,S,tminSST,beta,tmax,rmax);
toc;
%% Event Triggered Stabiliztation with Output Feedback
tic;
%S = ceil(2:m/8:m/2); % Sparsity levels
S = union(1,2:m/5:4*m/5);
tmax = 75; % 10sec of simulation time
beta = [0.5,0.5]; % beta_1, beta_2
x0 = 1*ones(n,1)/sqrt(n); % Initial state with l_2^2 2
% Compute Q2
Peff = P*F*C;
Q2 = 1*eigs(-(1/beta(1))*(Peff')*(Q0\Peff),1,'largestabs')*eye(n);
P2 = lyap(A1',-Q2);
tminSET = OutputFeedbackETStabilization(A,B,C,K,K_part,F,PO,P2,P,x0,S,beta,tmax);
toc;