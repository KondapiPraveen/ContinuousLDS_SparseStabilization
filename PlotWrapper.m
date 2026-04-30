% Plot Wrapper script
% CT-LDS : N mass-spring system 
% Find out more at http://www.ece.umn.edu/~mihailo/software/lqrsp/mass_spring.html

close all; clear; clc;
% N denotes number of masses in the mass-spring system
N = 4*10; c = zeros(N,1); n = 2*N;
% Entries in the main, sub/super diagonal in symmetric toeplitz matrix
c(1) = -2; c(2) = 1;
T = toeplitz(c);

% System matrices generation: A = [0 I;T 0], B = [0;I], C = [I 0] 
A = zeros(n,n); % Transfer matrix
A(1:N,N+1:n) = eye(N);
A(N+1:n,1:N) = T;

B = zeros(n,N); % Input matrix
% Effective Input Matrix
B_eff = eye(N);
%B_eff = diag(0.5*(1 + rand(N,1))); % Unif(0.5,1)
B(N+1:n,:) = B_eff; m = N;

p = N; C = zeros(p,n);
C(:,1:N) = eye(N);

% run PBH test to check if (A,B) is stabilizable
% Ths test is eecommended to run for new systems
checkStable = 1;
if checkStable
    isStable = isStabilizable(A,B);
    if isStable == 0
        error("The system is not stabilizable. Please choose a stabilizable system");
    end

    isDetect = isStabilizable(A.',C.');
    if isDetect
        error("The system is not detectable. Please choose a detectable system");
    end
end

% Compute the stabilizing (non-sparse) controller using lqr MATLAB function
% Find K s.t. A-BK is hurwitz using lqr function. Control input: u=-Kx
Q = eye(n); R = eye(m); % state and input regularization for state and input
K = lqr(A,B,Q,R); % Find the feedback matrix
F = lqr(A.',C.',Q,R); % Find the Luenberger observer gain

% Compute the lyapunov function V = x^{\T}Px s.t. (A-BK)^{\T}P + P(A-BK)
% We compute this using lyap fuction in MATLAB
K = -K; % In our work we consider A+BK is Hurwitz.
A0 = A+B*K; % Closed loop dynamics
Q0 = eye(n); % A0^{\T}*P + P*A0 + Q0 = 0
P = lyap(A0.',Q0); % lyap solves for X in A*X + X*A^{\T} + Q = 0, for A,Q inputs

% Compute Optimal Partition
K_row_l2 = vecnorm(K,2,2).^2; % Compute l2 norm of each row of K
[~,K_part] = sort(K_row_l2,'descend');
meff = m - sum(K_row_l2 == 0); % number of all zero rows in K
K_part_eff = K_part(1:meff);
%% State Feedback Continuous Control
% Initialization
tmax = 10; % 10sec of simulation time
S = 1:m/5:3*m/5;
x0 = ones(n,1)/sqrt(n/2); % Initial state with l_2^2 2
c4 = eigs(Q0,1,'smallestabs');
assert(c4 > 0,'System is not stabilizable via switching') % Otherwise system is not stabilizable via switching
gamma = c4/2; % gamma < c4
assert(c4 > gamma,'Non Zeno behavior is not guaranteed'); % To avoid Zeno behavior
% -- add a brace to skip this simulation
StateFeedbackContinuousControl(A,B,K,P,x0,K_part_eff,S,gamma,tmax);
%}
%% Feedback Gain under different partitions and convex combinations
% Choose Random B_eff
S = 1:m/5:m;
CtrlEnergyBound(K,K_part_eff,S);
%% Periodic Switching for Dynamic without and with Bounded Noise
S = 1:m/5:3*m/5;
rhoS = PeriodicSwitching(A,B,C,K,F,x0,K_part,S);
%% Periodic Switching for Static Output Feedback w/0, w/ Bounded Noise