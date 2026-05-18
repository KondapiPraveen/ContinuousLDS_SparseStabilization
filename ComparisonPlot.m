% Comparison of the our solution with the method from the papers: 
% Sparse control of linear continuous-time systems: A geometric approach

% System Matrices
clear; clc; close all;
gamma1 = 1;
A = [0,0,0;0,0,0;1,1,gamma1]; B = [1,0;0,1;0,0];
n = size(A,1); m = size(B,2); L = m;
x0 = ones(n,1)/sqrt(n); % Initial state vector
tmax = 20; % 10sec of simulation time

% Feedback matrices
k11 = 2;k13 = 1;k22 = 1;k23 = 1;
K1 = [(-k11-gamma1), 0, (-gamma1^2-gamma1*k11-k13); zeros(1,n)];
K2 = [zeros(1,n); 0, -k22-gamma1, -gamma1^2-gamma1*k22-k23];
Kl = zeros(m,n,L); Kl(:,:,1) = K1; Kl(:,:,2) = K2;
alpha = 0.5;

A1 = A+B*K1; A2 = A+B*K2; % Different Subsystems
Al = zeros(n,n,L); Al(:,:,1) = A1; Al(:,:,2) = A2;
Aalpha = alpha*A1 + (1-alpha)*A2;

Q0 = eye(n); c4 = eigs(Q0,1,'smallestabs'); gamma = c4/2;
P = lyap(Aalpha',Q0); % V(x) = x'*P*x

% Compute Q1, Q2 and the Lower Bound
Ql = zeros(n,n,L); Ql(:,:,1) = A1'*P + P*A1; Ql(:,:,2) = A2'*P + P*A2;
normAl = max(norm(A1),norm(A2));
nu1 = norm(A1'*(Ql(1)-Q0) + (Ql(1)-Q0)*A1); nu2 = norm(A2'*(Ql(2)-Q0) + (Ql(2)-Q0)*A2);
nu = max(nu1,nu2);
varTheta = 1 + (sqrt(1+8*((c4-gamma)*normAl/nu))-1)/4;
lwrBnd = min([(c4-gamma)/(varTheta^2*nu),log(varTheta)/normAl]);
Ql(:,:,1) = Ql(:,:,1) + gamma*eye(n); Ql(:,:,2) = Ql(:,:,2) + gamma*eye(n);

% Compute the trajectory and the control energy at each time instant using
% trajectory based switching
t = 0; x_curr = x0; cntr = 1; % Counter for storing x(t), t
maxEvents = ceil(tmax/lwrBnd)+1; % Max number of events for a given s
xl = zeros(maxEvents,1); t1 = zeros(maxEvents,1); ul = zeros(maxEvents,1);
xl(cntr) = (x_curr.')*x_curr; tl(cntr) = t; % mul with P for V(x)
cntr = cntr + 1;
Adelta = zeros(n,n,L);
for l = 1:L
    Adelta(:,:,l) = expm(Al(:,:,l)*lwrBnd);
end
while t <= tmax
    % Step 1: Find the subsystem
    l_curr = 1; l_curr_cost = (x_curr.')*Ql(:,:,1)*x_curr;
    for l = 2:L
        l_cost = (x_curr.')*Ql(:,:,l)*x_curr;
        if l_cost <= l_curr_cost
            l_curr = l;
            l_curr_cost = l_cost;
        end
    end

    % Compute the input energy
    ul(cntr-1) = norm(Kl(:,:,l_curr)*x_curr)^2;

    % Step 2: Find the next swtiching instant t_{k+1}
    t_next = t + lwrBnd;
    canItr = 1; xk = x_curr;
    while canItr
        xk = Adelta(:,:,l_curr)*xk;
        fval = (xk.')*Ql(:,:,l_curr)*xk;
        t_next = t_next + lwrBnd;
        if fval > 0
            canItr = 0;
            t_next = t_next - lwrBnd;
        end
    end
    assert(t_next > t, 't_next should be greater than t');
            
    % Step 3: Compute the state at t_{k+1} and update t, x_curr
    x_curr = expm(Al(:,:,l_curr)*(t_next-t))*x_curr;
    t = t_next;
    xl(cntr) = (x_curr.')*x_curr; tl(cntr) = t; cntr = cntr + 1;
end
%% Initialization for our method
% Compute the stabilizing (non-sparse) controller using lqr MATLAB function
% Find K s.t. A-BK is hurwitz using lqr function. Control input: u=-Kx
%
Q = eye(n); R = 6*eye(m); % state and input, cross regularization for state and input
%Q = eye(n); R = 2*eye(m); % state and input, cross regularization for state and input
N = zeros(n,m); % cross regularization term PS : DO NOT CHANGE
[K, P1] = lqr(A,B,Q,R,N); % Find the feedback matrix

% Compute the lyapunov function V = x^{\T}Px s.t. (A-BK)^{\T}P + P(A-BK)
% We compute this using lyap fuction in MATLAB
K = -K; % In our work we consider A+BK is Hurwitz.
%}
%K = alpha*K1 + (1-alpha)*K2;
A0 = A+B*K; % Closed loop dynamics
Q0 = eye(n); % A0^{\T}*P + P*A0 + Q0 = 0
P0 = lyap(A0',Q0); % lyap solves for X in A*X + X*A^{\T} + Q = 0, for A,Q inputs
%}

% Compute Optimal Partition
K_row_l2 = vecnorm(K,2,2).^2; % Compute l2 norm of each row of K
[~,K_part] = sort(K_row_l2,'descend');
meff = m - sum(K_row_l2 == 0); % number of all zero rows in K
K_part_eff = K_part(1:meff);
% State Feedback Continuous Control
% Initialization
tic;
S = 1;
c4 = eigs(Q0,1,'smallestabs');
assert(c4 > 0,'System is not stabilizable via switching') % Otherwise system is not stabilizable via switching
gamma = c4/2; % gamma < c4
assert(c4 > gamma,'Non Zeno behavior is not guaranteed'); % To avoid Zeno behavior
% Comment to skip this simulation
[tminS,Koptl,uoptl,toptl,xoptl,counter] = StateFeedbackContinuousControl(A,B,K,P0,x0,K_part_eff,S,gamma,tmax);
toc;
%% Plotting
figure
yyaxis left
str = "$s = 1$";
semilogy(toptl{1}(1:counter-1),xoptl{1}(1:counter-1),"DisplayName","Ours",'LineWidth',2.5,'Marker','s','MarkerSize',10);
xlabel('$time (t)$','Interpreter','latex','FontSize',16); ylabel('$\Vert \bf{x}(t) \Vert^2$','Interpreter','latex','FontSize',16);
%legend('interpreter','latex','FontSize',16,'NumColumns',2);
% plot for energy
yyaxis right
%figure();
semilogy(toptl{1}(1:counter-1),uoptl(1:counter-1),"DisplayName","Ours",'LineWidth',2.5,'Marker','d','MarkerSize',10); hold on
yyaxis left
semilogy(tl(1:cntr-1),xl(1:cntr-1),'--',"DisplayName","CDC'25",'LineWidth',2.5,'Marker','s','MarkerSize',10);
yyaxis right
semilogy(tl(1:cntr-1),ul(1:cntr-1),'--',"DisplayName","CDC'25",'LineWidth',2.5,'Marker','d','MarkerSize',10);
ylabel('$\Vert \bf{u}(t) \Vert^2$','Interpreter','latex','FontSize',16); %xlabel('$time (t)$','Interpreter','latex','FontSize',16);
legend('interpreter','latex','FontSize',16,'NumColumns',2); grid on;
%{
% Max frobenious norm
K_fro = max(norm(Kl(:,:,1),'fro'),norm(Kl(:,:,2),'fro'))
K_opt_fro = max(norm(Koptl{1}(:,:,1),'fro'),norm(Koptl{1}(:,:,2),'fro'))
%}