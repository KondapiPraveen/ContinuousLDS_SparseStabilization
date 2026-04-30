% Demonstrate energy usage for various s and various paritions and cvx comb
%{
Inputs:
    K  : Non-sparse stabilizing control, i.e., A0 = A+BK is Hurwitz
    K_part : Ordering of rows of K to obtain K_l's
    S : set of sparsity levels
%}
%{
    Plot 2: \max_l |K_l|_F^2 vs 's' for diff \alpha, \set{S_l}.
%}

function CtrlEnergyBound(K,K_part,S)
    % Initialization
    lenS = numel(S); meff = numel(K_part);
    EngBnd = zeros(4,lenS);
    UppBnd = zeros(1,lenS);

    for i=1:lenS
        s = S(i); L = ceil(meff/s); % Number of controllers for a given s
        K_row_l2 = vecnorm(K,2,2).^2;

        % S1: Optimal cvx comb. and partition
        CONST = 0; % Find the max. frobenious norm of feedback matrices
        for l=1:L
            Sl = K_part((l-1)*s+1:min(l*s,meff));
            CONST = CONST + sqrt(sum(K_row_l2(Sl)));
            if l==1
                UppBnd(1,i) = (L*CONST)^2;
            end
        end
        EngBnd(1,i) = CONST^2;
        
        % S2: Trivial Partition and Optimal cvx comb for the trv partition
        CONST=0; % Find max. frobenious norm of feedback matrices of given parition
        for l=1:L
            Sl = (l-1)*s+1:min(l*s,meff);
            CONST = CONST + sqrt(sum(K_row_l2(Sl)));
        end
        EngBnd(2,i) = CONST^2;

        % S3: Optimal partitiion and uniform cvx combination
        Sl = K_part(1:min(1*s,meff));
        EngBnd(3,i) = L^2*sum(K_row_l2(Sl));

        % S4: Uniform cvx Combination and Trivial Partition
        CONST = 0;
        for l=1:L
            Sl = (l-1)*s+1:min(l*s,meff);
            FeedFroNorm = sum(K_row_l2(Sl));
            if CONST < FeedFroNorm
                CONST = FeedFroNorm;
            end
        end
        EngBnd(4,i) = L^2*CONST;
    end

    % Plotting
    figure();
    Labels=["$\bs{alpha}^*, \mathcal{S}^*$","$\bs{alpha}^*, \mathcal{S}$","$\bs{alpha}, \mathcal{S}^*$","$\bs{alpha}, \mathcal{S}$"];
    Markers = {'s','d','o','+','x','*','^','v'}; % For plotting
    for i=1:lenS
        plot(S,EngBnd(i,:),'DisplayName',Labels(i),'Marker',Markers(i),'MarkerSize',10,'LineWidth',2.5); hold on;
    end
    plot(S,UppBnd,'DisplayName','Upper bound for $\bs{alpha}^*, \mathcal{S}^*$','ko-.','MarkerSize',10,'LineWidth',2.5); hold off
    legend('Interpreter','latex','NumColumns',2);
    ylabel("Sparsity ($s$)",'Interpreter','latex')
    xlabel("Max Sq. Frob. Norm $\max_l ||K_l||_F^2$",'Interpreter','latex');
end