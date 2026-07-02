%% ECC 2024: Stabilization of Diagonal Matrices - Script File
% Result 1 with diagonal system matrices
% Continuous feedback with precompuation of next switching time
%{
Plot 1 : ||x(t)|| over time (t) for various s
Plot 2 : ||u(t)|| over time (t) for various s
Plot 2b : ||K_l||/||K|| over time (t) for various s
Plot 3 : Support(u(t)) over time (t) for s=1
%}

%{
System Model
dot{x} = Ax + Bu
%}

clear; clc; close all; rng(0);
% Zacharay Karate Club Dataset import
load('./zachary.mat');

% System Matrices
% Use either zachary_bin or zachary_val as an argument for eig
[~,A] = eig(zachary_bin); n = size(A,2);
B = eye(n); % Input matrix
A = A/norm(A);

% Problem Parameters
S = [1,3,5,7]; lnS = length(S); % Set of sparsity levels
kappa = 0.06; gam = 0.02; % Specified rate of stability and threshold
t_init = 0; t_final = 25; % 8sec of Simulation
num_alpha = 2; % 1st is Optimal, 2nd is Uniform

% Initialization
x_0 = zeros(n,1); %x_0 = x_0/norm(x_0);
%x_0 = (1:n)'; x_0 = x_0/norm(x_0);
normalized_feedback = ones(num_alpha,lnS+1);
x_traj = cell(lnS+1,num_alpha); % Store the norm of the state trajectories
u_traj = cell(lnS+1,num_alpha); % Store the norm of the inputs
pre_time = t_init:0.1:t_final;
x_des = exp(-gam*pre_time)*1;
%x_des = exp(-kappa*pre_time)*1;
time = cell(lnS+1,num_alpha);
Supp_u = cell(num_alpha); % store the support for s=1


% Compute the entries to be controlled
sat_entries = kappa+diag(A);
N = 1:n; % Enumeration of first 'n' natural numbers 
Idx = N(sat_entries >= kappa - gam/2); % Indices to be controlled
len_Idx = numel(Idx); % Full Sparsity level
S = union(S,len_Idx);

x_0(Idx) = randn(numel(Idx),1); x_0 = x_0/norm(x_0); 

% Compute the all gains of the controlled indices
bar_b = -sat_entries(Idx)./diag(B(Idx,Idx));

% Compute the Feedback matrix for unconstrained problem
K_full = zeros(n);
K_full(Idx,Idx) = diag(bar_b); K_full_norm = norm(K_full);

for j=1:num_alpha
    for i=1:lnS
        s = S(i);
        % Design of switching subsystems
        L = ceil(numel(Idx)/s); % Number of subsystems
        A_l = zeros(n,n,L); K_l = zeros(n,n,L);
        
        % Compute values for optimal alpha
        abs_bar_b = zeros(L,1);
        for l=1:L
            start_idx = max(1,len_Idx-l*s+1);
            abs_bar_b(l) = max(abs(bar_b(start_idx:end-(l-1)*s)));
        end
        sum_abs_bar_b = sum(abs_bar_b,"all");

        % Compute alpha
        if j==1
            cvx_comb = abs_bar_b/sum_abs_bar_b; % Optimal alpha
        else
            cvx_comb = ones(L,1)/L; % Uniform alpha
        end
        
        % Compute A_l = A+B*K_l
        max_K_l = 0;
        for l=1:L
            u = zeros(n,1);
            start_idx = max(1,len_Idx-l*s+1);
            u(Idx(start_idx:end-(l-1)*s)) = bar_b(start_idx:end-(l-1)*s)/cvx_comb(l);
            K_l(:,:,l) = diag(u);
            A_l(:,:,l) = A+B*K_l(:,:,l);
            curr_K_l = max(abs(u));
            if max_K_l < curr_K_l
                max_K_l = curr_K_l;
            end
        end
        normalized_feedback(j,i) = norm(max_K_l)/K_full_norm;

        % Start of Simulation of Continuous Dynamical systems
        t=0; r = 1; % time and # switchings
        time{i,j}(r) = 0; 
        x_curr = x_0;
        x_traj{i,j}(r) = norm(x_curr);
        while t<t_final
            % Find the correct A_l
            min_sub = 1; % Initialization for correct subsystem
            Vx = zeros(n,1); % Compute the Lyapunov function derivative
            for l=1:L
                Vx(l) = (x_curr.')*A_l(:,:,l)*x_curr;
                if l > 1 && Vx(l) < Vx(min_sub)
                    min_sub = l;
                end
            end
            % Compute the control energy
            u_traj{i,j}(r) = norm(K_l(:,:,min_sub)*x_curr);
            %u_traj{i}(r) = norm(K_l(:,:,min_sub));
    
            % Store the support
            if s==1
                Supp_u{j}(r) = Idx(L-min_sub+1);
            end
    
            % Find t_r using newton raphson
            canItr = 1; t_next = t;
            while canItr
                t_prev = t_next;
                v = (x_curr.')*diag(exp((diag(A_l(:,:,min_sub))+(gam/2))*(t_next - t)));
                t_next = t_next - ((v*(2*A_l(:,:,min_sub) + gam*eye(n))*(v.'))/(v*((2*A_l(:,:,min_sub) + gam*eye(n)).^2)*(v.')));
                canItr = (abs(t_next-t_prev)/abs(t_prev) > 1e-3);
            end
            r = r+1;
            x_curr = exp(diag(A_l(:,:,min_sub))*(t_next - t)).*x_curr;
            x_traj{i,j}(r) = norm(x_curr).^2;
            time{i,j}(r) = t_next; t = t_next;
        end
        %}
    end
end

% Compute 
time{lnS+1} = t_init:0.5:t_final + 0.5;
x_full = exp(diag(A + B*K_full)*(time{lnS+1})).*repmat(x_0,[1,length(time{lnS+1})]);
x_traj{lnS+1} = vecnorm(x_full).^2;
u_traj{lnS+1} = vecnorm(K_full*x_full(:,1:end-1)); clear x_full;

%% Plotting
DefaultColors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; 0.3010 0.7350 0.9330; 0.6350 0.0780 0.1840];

% Plot 1 : ||x(t)|| vs t for various 's'
figure();
markers = ['s','d','o','x','+','*'];
mk_idx = 1:10:length(pre_time);
semilogy(pre_time,x_des,'k--','Linewidth',2.5,'DisplayName','Upper Bound','MarkerSize',10,'Marker',markers(lnS+1),'MarkerIndices',mk_idx);
hold on;

set(gca(),'FontSize',16,'FontWeight','normal');
for i=lnS+1:-1:1
    str = strcat('$s$=',num2str(S(i)));
    mk_idx = 1:25:length(time{i,1});
    semilogy(time{i,1},x_traj{i,1},'Linewidth',2.5,'DisplayName',str,'MarkerSize',10,'Marker',markers(i),'MarkerIndices',mk_idx,'Color',DefaultColors(i,:));
    legend('Interpreter','latex','FontSize',14,'NumColumns',2);
    hold on;
    xlabel('Time $(t)$','Interpreter','latex');
    str = sprintf('State vector energy $(||x(t)||^2)$');
    ylabel(str,'Interpreter','latex');
    xlim([0,t_final])
end

%% Plot 2a : ||u(t)|| vs t for various 's'
figure();
for i=1:lnS+1
    mk_idx = 1:5:(length(time{i,1})-1);
    str = strcat('$s$=',num2str(S(i)));
    semilogy(time{i,1}(1:end-1),u_traj{i},'Linewidth',2.5,'DisplayName',str,'MarkerSize',10,'Marker',markers(i),'MarkerIndices',mk_idx,'Color',DefaultColors(i,:));
    hold on;
end
set(gca(),'FontSize',16,'FontWeight','normal');
legend('Interpreter','latex');
xlabel('Time $(t)$','Interpreter','latex');
str = sprintf('Input energy $(||u(t)||)$');
ylabel(str,'Interpreter','latex');
%%
% Plot 2b : ||K_l||/||K|| vs sparsity 's'
figure();
S_sup = union(S,len_Idx);
upp_bnd = ceil((len_Idx./S_sup));
loglog(S_sup,normalized_feedback(2,:),'LineWidth',2.5,'DisplayName','$\alpha = \frac{\textbf{1}}{L}$','Marker','s','MarkerSize',10,'Color',DefaultColors(1,:)); hold on
loglog(S_sup,upp_bnd,'k--o','LineWidth',2.5,'DisplayName','Upper Bound','MarkerSize',10);
loglog(S_sup,normalized_feedback(1,:),'LineWidth',2.5,'DisplayName','$\alpha = \alpha^*$','Marker','s','MarkerSize',10,'Color',DefaultColors(2,:));
set(gca(),'FontSize',16,'FontWeight','normal');
legend('Interpreter','latex');
xlabel('Sparsity $(s)$','Interpreter','latex')
xticks([1,3,5,7,15])
yticks([1,2,5,8,15])
str = 'Feedback Gain $(\max_l \hat{k}_l/k^*)$';
ylabel(str,'Interpreter','latex','FontSize',18);
%%
% Plot 3: Support of u(t)
fig = figure();
if S(1) == 1
    subplot(2,1,1)
    stairs(time{1}(1:end-1),Supp_u{1},'Linewidth',2,'MarkerSize',10);
    xlim([0,6])
    set(gca(),'FontSize',14,'FontWeight','normal');
    subplot(2,1,2)
    stairs(time{1}(1:end-1),Supp_u{1},'Linewidth',2,'MarkerSize',10);
    xlim([t_final-6,t_final])
    set(gca(),'FontSize',14,'FontWeight','normal');
    % Common labels
    ax = axes(fig,'visible','off');
    ax.XLabel.Visible = 'on';
    ax.YLabel.Visible = 'on';
    ylabel(ax,'Eigenvalue indices','FontSize',18,'Interpreter','latex')
    xlabel(ax,'Time $(t)$','Interpreter','latex','FontSize',18);
    %exportgraphics(fig,'R1_Plt3_USupp_Time.pdf')
end
