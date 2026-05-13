% State Feedback based Self Triggered Sparse Stabilization of LDS
%{
Inputs
    System Matrices : A,B
    Feedback Matrix : K
    Ordering of rows of K : K_part
    Lyapunov Matrix: P
    Inital state : x0
    Sparsity : S
    Inter Event Time lower bounds : tminSST
    Max. simulation time: tmax
    max time of non-interrupted control: rmax
%}
function STSparseStabilization(A,B,K,K_part,P,x0,S,tminSST,beta,tmax,rmax)
    lenS = numel(S); meff = numel(K_part); n = size(A,1);
    Fs = zeros(2*n,2*n); Fs(1:n,1:n) = beta*(A+B*K); % Auxilary system
    figure(); Markers = {'s','d','o','+','x','*','^','v'}; % For plotting
    legend('NumColumns',2,'Interpreter','latex');

    for i=1:lenS
        s = S(i); % Sparsity
        Del = 1*tminSST(i); % Lower bound on the inter-event for s
        MaxSwitches = ceil(tmax/Del);
        ExpFs = expm(Fs*Del);
        
        % Initialization
        L = ceil(meff/s); % # Subsystems
        Fl = zeros(2*n,2*n,L); % Dynamics of [x, x_e]
        ExpFl = zeros(2*n,2*n,L); % Matrix exponential of Fl's
        Ql = zeros(n,n,L); %  Matrix in the derivative of V(t)
        xNorm2 = zeros(MaxSwitches+1,1); xNorm2(1) = norm(x0);
        time = zeros(MaxSwitches+1,1);
        
        % Obtain optimal feedback matrices and cvx combinations
        [Kl,~] = OptimalFeedback(K,K_part,s); % Get Optimal Feedback and Cvx. Comb.

        % Compute the subsystems Fl's
        for l=1:L
            Beff = B*Kl(:,:,l); Aeff = A+Beff;
            Fl(:,:,l) = [Aeff, Beff; -Aeff, -Beff];
            ExpFl(:,:,l) = expm(Fl(:,:,l)*Del);
            % Compute Ql = (A+B*K_l)'*P + P*(A+B*K_l)
            Ql(:,:,l) = Aeff'*P + P*Aeff;
        end
        
        t = 0; % Current time
        count = 1; zeroVec = zeros(n,1);
        x_ext = [x0;zeroVec]; % [x,x_e]
        while t < tmax
            s_ext = x_ext; % The state after update is same of both sys.

            % Compute the appropriate subsystem
            l_curr = 1; l_curr_cost = (x_ext(1:n)')*Ql(:,:,1)*x_ext(1:n); 
            for l=1:L
                l_cost = (x_ext(1:n)')*Ql(:,:,l)*x_ext(1:n);
                if l_cost <= l_curr_cost
                    l_curr = l;
                    l_curr_cost = l_cost;
                end
            end
            
            x_ext = ExpFl(:,:,l_curr)*x_ext;
            s_ext = ExpFs*s_ext; t = t + Del;
            for k=1:rmax
                x_ext_next = ExpFl(:,:,l_curr)*x_ext;
                s_ext_next = ExpFs*s_ext;
                
                Vt = (x_ext_next(1:n)')*P*x_ext_next(1:n); % V(t)
                St = (s_ext_next(1:n)')*P*s_ext_next(1:n); % S(t)

                if Vt > St
                    break;
                end
                t = t + Del;
                x_ext = x_ext_next;
                s_ext = s_ext_next;
            end
            x_ext = [x_ext(1:n); zeroVec]; count = count + 1;
            xNorm2(count) = norm(x_ext(1:n));
            time(count) = t;
        end

        xNorm2 = xNorm2.^2; % Get the sq. l2 norm of the state
        % Plotting
        str = sprintf('$s = %d$',s);
        semilogy(time(1:count),xNorm2(1:count),'DisplayName',str,'Marker',Markers(i),'MarkerSize',10,'LineWidth',2.5);
        hold on;
    end
    legend('NumColumns',2,'Interpreter','latex','FontSize',16); 
    xlabel('time (s)','Interpreter','latex','FontSize',16); ylabel('$\Vert \textbf{x}(t) \Vert^2$','Interpreter','latex','FontSize',16);
    title('Self Triggered State Feedback')
end