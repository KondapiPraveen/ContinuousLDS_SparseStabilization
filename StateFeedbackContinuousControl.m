% State Feedback based Continuous Sparse Stabilizing Control
%{ 
    Problem:
    \dot{x}(t) = Ax(t) + Bu(t)
    u(t) = K_{\sigma(t)}x(t)
%}

%{
Inputs:
    A,B: System matrices (transfer, input)
    K  : Non-sparse stabilizing control, i.e., A0 = A+BK is Hurwitz
    P  : Lyapunov matrix, i.e., A0^{\T}P + P*A0 is ND
    x0 : initial state
    K_part : Ordering of rows of K to obtain K_l's
    S : set of sparsity levels
    gamma : rate of stability
    tmax : Max sim. time 0 <= t <= tmax
%}

%{
Outputs:
    Plot 1: |x(t)| vs t for various 's'
    Table 1: Lower min inter event time for various s
%}

function [tmin,Kl,ul,tl,xl,cntr] = StateFeedbackContinuousControl(A,B,K,P,x0,K_part,S,gamma,tmax)
    % Initialization
    Q = ((A+B*K).')*P + P*(A+B*K); c4 = eigs(-Q,1,'smallestabs');
    lenS = numel(S); n = size(A,1); m = size(B,2); meff = numel(K_part);
    tmin = zeros(lenS,1); matI = eye(n);
    Al = cell(lenS,1); Ql = cell(lenS,1); tQl = cell(lenS,1); %tQl2 = cell(lenS,1);
    xl = cell(lenS,1); Kl = cell(lenS,1); tl = cell(lenS,1);
    tx = linspace(0,tmax,20); c3 = norm(P);
    UppBnd = x0'*P*x0*exp(-gamma*tx/c3);
    for i=1:lenS
        s = S(i); L = ceil(meff/s);
        Al{i} = zeros(n,n,L); % A_l : Closed loop dynamics of each subsystem
        Ql{i} = zeros(n,n,L); % Q_l = A_l.'*P + P*A_l 
        tQl{i} = zeros(n,n,L); % \tilde{Q}_l = A_l.'*(Q_l + \gamma*I) + (Q_l + \gamma*I)*A_l.'
        %tQl2{i} = zeros(n,n,L); % \tilde{Q}_l = A_l.'*(Q_l - Q) + (Q_l - Q)*A_l.'
    end
    figure;
    Markers = {'s','d','o','+','x','*','^','v'}; % For plotting
    
    for i = 1:lenS
        s = S(i); L = ceil(meff/s); % Number of controllers for a given s
        K_row_l2 = vecnorm(K,2,2).^2;
        Adelta = zeros(n,n,L);
        CONST = 0; % Find the optimal combination constant
        for l=1:L
            Sl = K_part((l-1)*s+1:min(l*s,meff));
            CONST = CONST + sqrt(sum(K_row_l2(Sl)));
        end
        
        % Design of A_l's where Al = A + B*Kl, Q_l, \tilde{Q}_l
        lwrBnd1 = 0; lwrBnd2 = 0;
        for l = 1:L
            Sl = K_part((l-1)*s+1:min(l*s,meff));
            K_eff = zeros(m,n);
            K_eff(Sl,:) = K(Sl,:);
            %tol = max(abs(K_eff),[],"all");
            %K_eff(abs(K_eff) < tol*1e-3) = 0;
            K_eff = (CONST/norm(K_eff,'fro'))*K_eff;
            Kl{i}(:,:,l) = K_eff;
            Al{i}(:,:,l) = A + B*K_eff;
            Ql{i}(:,:,l) = (Al{i}(:,:,l).')*P + P*Al{i}(:,:,l);
            tQl{i}(:,:,l) = (Al{i}(:,:,l).')*(Ql{i}(:,:,l)+gamma*eye(n)) + (Ql{i}(:,:,l)+gamma*eye(n))*(Al{i}(:,:,l).');
            %tQl2{i}(:,:,l) = (Al{i}(:,:,l).')*(Ql{i}(:,:,l)-Q) + (Ql{i}(:,:,l)-Q)*(Al{i}(:,:,l).');
            normAl = sqrt(eigs((Al{i}(:,:,l).')*Al{i}(:,:,l),1));
            nu = abs(eigs((Al{i}(:,:,l).')*(Ql{i}(:,:,l)-Q) + (Ql{i}(:,:,l)-Q)*(Al{i}(:,:,l).'),1));
            if normAl > lwrBnd1
                lwrBnd1 = normAl;
            end
            if nu > lwrBnd2
                lwrBnd2 = nu;
            end
            Ql{i}(:,:,l) = Ql{i}(:,:,l) + gamma*eye(n);
        end

        % Compute the lower bound
        varTheta = 1 + (sqrt(1+8*((c4-gamma)*lwrBnd1/lwrBnd2))-1)/4;
        lwrBnd = min([(c4-gamma)/(varTheta^2*lwrBnd2),log(varTheta)/lwrBnd1]);
        tmin(i) = lwrBnd; % Save the lower bound for TABLE
        maxEvents = ceil(tmax/lwrBnd)+1; % Max number of events for a given s
        xl{i} = zeros(maxEvents,1); tl{i} = zeros(maxEvents,1);
        ul = zeros(maxEvents,1);
        for l = 1:L
            Adelta(:,:,l) = expm(Al{i}(:,:,l)*lwrBnd);
        end

        % Switching Stabilization
        t = 0; x_curr = x0; cntr = 1; % Counter for storing x(t), t
        xl{i}(cntr) = (x_curr.')*P*x_curr; % Change to P for V(x)
        tl{i}(cntr) = t; cntr = cntr + 1;
        while t <= tmax
            % Step 1: Find the subsystem
            l_curr = 1; l_curr_cost = (x_curr.')*Ql{i}(:,:,1)*x_curr;
            for l = 2:L
                l_cost = (x_curr.')*Ql{i}(:,:,l)*x_curr;
                if l_cost <= l_curr_cost
                    l_curr = l;
                    l_curr_cost = l_cost;
                end
            end

            % Compute the input energy
            ul(cntr-1) = norm(Kl{i}(:,:,l_curr)*x_curr)^2;

            % Step 2: Find the next swtiching instant t_{k+1}
            t_next = t + lwrBnd;
            canItr = 1; xk = x_curr;
            while canItr
                xk = Adelta(:,:,l_curr)*xk;
                fval = (xk.')*Ql{i}(:,:,l_curr)*xk;
                t_next = t_next + lwrBnd;
                if fval > 0
                    canItr = 0;
                    t_next = t_next - lwrBnd;
                end
            end
            %{
            while canItr
                xk = expm(Al{i}(:,:,l_curr)*(t_next-t))*x_curr;
                fval = (xk.')*Ql{i}(:,:,l_curr)*xk;
                if fval < 0
                    delta = 2*delta; t_next = t + delta;
                else
                    t1 = t + delta/2;
                    t2 = t_next;
                    canItr = 0;
                end
            end
            canItr = 1;
            while canItr
                t_next = (t1 + t2)/2;
                xk = expm(Al{i}(:,:,l_curr)*(t_next-t))*x_curr;
                fval = (xk.')*Ql{i}(:,:,l_curr)*xk;
                if fval < 0
                    t1 = t_next;
                else
                    t2 = t_next;
                end
                if t2 - t1 < lwrBnd
                    canItr = 0;
                end
            end
            %}
            assert(t_next > t, 't_next should be greater than t');
            
            % Step 3: Compute the state at t_{k+1} and update t, x_curr
            x_curr = expm(Al{i}(:,:,l_curr)*(t_next-t))*x_curr;
            t = t_next;
            xl{i}(cntr) = (x_curr.')*P*x_curr; % Change to P for V(x)
            tl{i}(cntr) = t; cntr = cntr + 1;
        end
        % Plotting : |x(t)| vs t for various 's'
        semilogy(tl{i}(1:cntr-1),xl{i}(1:cntr-1),'LineWidth',2.5,'Marker',Markers{i},'MarkerSize',10);
        hold on;
    end
    semilogy(tx,UppBnd,'k--o','LineWidth',2.5,'MarkerSize',10,'MarkerIndices',1:2:20); hold on
    xlabel('time $(t)$','Interpreter','latex','FontSize',16); ylabel('$V(t)$','Interpreter','latex','FontSize',16);
    legend([strcat('$s = $',num2str(S.'));"Upper Bound"],'interpreter','latex','FontSize',16,'NumColumns',2); grid on;
end