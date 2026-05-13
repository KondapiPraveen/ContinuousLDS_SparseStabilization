% Compute the SOF matrix for the LDS
%{
Inputs:
    System Matrices: A,B,C
    LQR cost matrices: Q, R, N
    Lyapunov matrix: P,
    maximum nubmer of Iterations: maxItr (choose sufficiently)
    Threshold for stopping criteria: thrsh (choose a small value)
%}
%{
Output: SOF matrix: F
%}

function F = SOF(A,B,C,Q,R,N,P,maxItr,thrsh)
    itr = 0;
    critVal = inf;
    pinvC = pinv(C);
    while (itr<maxItr) && (critVal > thrsh)
        K = R\(B'*P+N');
        F = K*pinvC;
        Atilde = A-B*F*C;
        Ftilde = F*C;
        Ntilde = N*Ftilde;
        RE = Atilde'*P + P*Atilde + Q + Ftilde'*R*Ftilde - Ntilde' - Ntilde;
        P = P + lyap((A-B*K)',RE);
        critVal = trace(RE^2);
        itr = itr + 1;
    end
end