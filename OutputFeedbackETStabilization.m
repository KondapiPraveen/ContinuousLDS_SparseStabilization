% Compute the Lower Bound on the Inter Event Times for output feedback
% based Event Triggered Dynamic stabilization

%{
Inputs:
    System Matrices: A,B,C
    Feedback Matrix: K
    Row Ordering of K: K_part
    Luenberger Gain: F
    Lyapunov Matrices: P1, P2, P such that
    (A+FC)'*P1 + P1*(A+FC) = -I
    (A+FC)'*P2 + P2*(A+FC) = Q2
    (A+BK)'*P + P*(A+BK) = Q
    initial state: x0
    Sparsity: S
    beta: [beta_1,\beta_2]
    Max. Sim. time: tmax
%}

% Output: Lower Bound on Inter-Event Times

function tmin = OutputFeedbackETStabilization(A,B,C,K,K_part,F,P1,P2,P,x0,S,beta,tmax)
    % Precomputation
    Ceff = F*C; ACeff = A+Ceff; xn = norm(x0); n = size(A,1);
    Q2 = (ACeff.')*P2 + P2*ACeff; Q = ((A+B*K).')*P + P*(A+B*K);
    Q3 = [beta(1)*Q P*Ceff; (P*Ceff).' Q2]; nQ3 = eigs(Q3,1,'largestabs');
    nQ = eigs(Q,1,'largestabs'); P1max = norm(P1); P1min = eigs(P1,1,'smallestabs'); P1maxmin = sqrt(P1max/P1min);
    matI = eye(2*n); tx = linspace(0,tmax,20);
    %Del = 1e-3;

    % Initialization
    lenS = numel(S); n = size(A,1); meff = numel(K_part); tmin = zeros(lenS,1);
    P3 = zeros(2*n,2*n);
    P3(1:n,1:n) = P; P3(n+1:2*n,n+1:2*n) = P2; S1 = [eye(n) zeros(n)];
    P3min = eigs(P3,1,'smallestabs'); P3max = eigs(P3,1,'largestabs');
    Aext = [A zeros(n) Ceff; -A zeros(n) -Ceff; zeros(n) zeros(n) ACeff];


    fig1 = figure; %fig2 = figure;
    Markers = {'s','d','o','+','x','*','^','v'}; % For plotting
    h1 = @(t) P1maxmin*xn*exp(-t/(2*P1max));
    
    c10 = P3max*P1maxmin; % Divide by P3min for ||Psi(t)||^2
    c9 = (1-beta(2))*nQ3/P3max;
    c11 = max(c9,-1/P1max);
    g1 = @(t) 2*c10*xn^2*exp(c11*t);
    UppBnd = max(2*P3max*h1(tx).^2,g1(tx)); % Remove P3max for ||Psi(t)||^2

    for i=1:lenS
        s = S(i); L = ceil(meff/s); % Spasity and the # Subsystems

        % Initialization
        Fl = zeros(3*n,3*n,L); % Dynamics of [x_hat, x_hat_e, e]
        ExpFl = Fl; % Store matrix exponentials of Fl
        Al = zeros(2*n,2*n,L); % Transfer Matrix of [x_hat, e]
        Ql = zeros(n,n,L); Qs = zeros(n,n,L);
        Bl = zeros(2*n,n,L); % Effective Weight for Deterministic Noise (x_hat_e)
        
        % Optimal Feedback matrices
        [Kl,~] = OptimalFeedback(K,K_part,s); % Get Optimal Feedback and Cvx. Comb.

        a = 0; b = 0; c1 = 0;
        for l=1:L
            Beff = B*Kl(:,:,l); Aeff = A+Beff;
            Fl(:,:,l) = [Aeff Beff Ceff; -Aeff -Beff -Ceff; zeros(n) zeros(n) ACeff];
            Al(:,:,l) = [Aeff Ceff; zeros(n) ACeff];
            Bl(:,:,l) = [Beff; zeros(n)];
            Ql(:,:,l) = Aeff'*P + P*Aeff;
            Qs(:,:,l) = Ql(:,:,l) - beta(1)*Q;

            al = norm(Al(:,:,l)); bl = norm(Bl(:,:,l)); cl = norm(P*Beff);
            a = al*(al > a) + a*(al <= a);
            b = bl*(bl > b) + b*(bl <= b);
            c1 = cl*(cl > c1) + c1*(cl <= c1);
        end

        c7 = -beta(2)*nQ3/(2*c1); c2=0; 
        for l=1:L
            Qtil = S1'*Qs(:,:,l)*S1;
            Atil = Qtil*Al(:,:,l);
            cl = norm(Atil + Atil') + 2*norm(Qtil*Bl(:,:,l))*c7;
            c2 = cl*(cl > c2) + c2*(cl <= c2);
        end

        c8 = 2*P3max*c2/P3min;
        % Lower Bound on Inter-Event Time
        LowerBound1 = -(1-beta(1))*nQ/c8;
        LowerBound2 = log((a*(1+c7))/(a+b*c7))/(a-b);
        LwrBnd = min(LowerBound2,LowerBound1);
        tmin(i) = LwrBnd;

        MaxSwitches = ceil(tmax/LwrBnd);
        xNorm2 = zeros(MaxSwitches+2,1); %xNorm2(1) = 0; 
        xNorm2(1) = [zeros(n,1); -x0]'*P3*[zeros(n,1); -x0]; % change matI to P3 to plot V(t)
        Enorm2 = zeros(MaxSwitches+2,1); Enorm2(1) = xn;
        time = zeros(MaxSwitches+2,1);
        ExpAext = expm(Aext*LwrBnd);
        %ExpAext = expm(Aext*Del);
        
        for l=1:L
            ExpFl(:,:,l) = expm(Fl(:,:,l)*LwrBnd);
            %ExpFl(:,:,l) = expm(Fl(:,:,l)*Del);
        end

        % Start Simulation: Check condtions at intervals of LwrBnd
        t=0; % Current time
        canUpdate = 1; % Can I update the input
        count = 1;
        x_ext = [zeros(n,1);zeros(n,1);-x0];
        while canUpdate
            x_ext = ExpAext*x_ext; % state at next period
            x_hat = x_ext(1:n); % Extract \hat{x}
            if norm(x_hat) >= h1(t)
                canUpdate = 0;
                % Compute l_curr
                l_curr = 1; l_curr_cost = Inf; % Current active subsystem
                for l=1:L
                    l_cost = x_hat'*Ql(:,:,l)*x_hat;
                    if (l_cost < l_curr_cost)
                        l_curr = l;
                        l_curr_cost = l_cost;
                    end
                end
                % Update the input
                x_ext(n+1:2*n) = zeros(n,1);
                % Store x(t)
                x_curr = [x_hat; x_ext(2*n+1:3*n)];
                %x_curr = x_hat;
                count = count + 1;
                %xNorm2(count) = norm(x_curr);
                xNorm2(count) = x_curr'*P3*x_curr; % change matI to P3 to plot V(t)
                Enorm2(count) = norm(x_ext(2*n+1:3*n));
                time(count) = t + LwrBnd;
            end
            % Update time
            t = t + LwrBnd;
            %t = t + Del;
        end

        while t < tmax
            x_ext = ExpFl(:,:,l_curr)*x_ext;
            x_hat = x_ext(1:n); x_hat_e = x_ext(n+1:2*n);
            
            % Check for update conditions
            quant1 = x_hat'*Qs(:,:,l_curr)*x_hat;
            quant2 = norm(x_hat_e)/norm(x_hat);
            if (quant1 > 0 || quant2 > c7) && norm(x_hat) >= h1(t)
                % Compute l_curr
                l_curr = 1; l_curr_cost = Inf; % Current active subsystem
                for l=1:L
                    l_cost = x_hat'*Ql(:,:,l)*x_hat;
                    if (l_cost < l_curr_cost)
                        l_curr = l;
                        l_curr_cost = l_cost;
                    end
                end
                % Update the input
                x_ext(n+1:2*n) = zeros(n,1);
                % Store x(t)
                x_curr = [x_hat; x_ext(2*n+1:3*n)];
                %x_curr = x_hat;
                count = count + 1;
                %xNorm2(count) = norm(x_curr);
                Enorm2(count) = norm(x_ext(2*n+1:3*n));
                xNorm2(count) = x_curr'*P3*x_curr; % change matI to P3 to plot V(t)
                time(count) = t + LwrBnd;
            end
            % Update time
            t = t + LwrBnd;
            %t = t + Del;
        end

        % Plot
        %xNorm2 = xNorm2.^2;
        figure(fig1)
        str = sprintf('$s = %d$',s);
        MarkerInd = floor(linspace(1,count,10));
        semilogy(time(1:count),xNorm2(1:count),'DisplayName',str,'Marker',Markers(i),'MarkerIndices',MarkerInd,'MarkerSize',10,'LineWidth',2.5);
        hold on;
        
        %{
        Plot for debugging
        figure(fig2)
        Enorm2 = Enorm2.^2;
        semilogy(time(1:count),Enorm2(1:count),'DisplayName','Error','LineWidth',2.5); hold on
        semilogy(time(1:count),h1(time(1:count)).^2,'DisplayName','h(t)','LineWidth',2.5); hold off
        legend('NumColumns',2,'Interpreter','latex','FontSize',16);
        %}
    end
    figure(fig1)
    semilogy(tx,UppBnd,'k--o','DisplayName',"Upper Bound",'MarkerSize',10,'LineWidth',2.5); grid on; hold off
    legend('NumColumns',2,'Interpreter','latex','FontSize',16); 
    xlabel('time $(t)$','Interpreter','latex','FontSize',16); %ylabel('$\Vert \textbf{x}(t) \Vert^2$','Interpreter','latex','FontSize',16);
    %ylabel('$\Vert \psi(t) \Vert^2$','Interpreter','latex','FontSize',16);
    ylabel('$V(t)$','Interpreter','latex','FontSize',16);
    title('Noiseless Dynamic Feedback')
end