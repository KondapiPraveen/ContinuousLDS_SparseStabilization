% Function to test whether (A,B) is stabilizable
%{
    Input : System matrices A, B
    Ouput : isStable \in {0,1}
    isStable = 1 if the system is stabilizable otherwise 0
%}
function isStable = isStabilizable(A,B)
    n = size(A,1); m = size(B,2);
    eigVals = eig(A); % All of the eigenvalues of the matrix A
    impEigVals = eigVals(real(eigVals) >= 0); % Candidate eigenvalues for PBH test
    isStable = 1; % Variable to identify stabilizability
    for eVal = impEigVals.'
        testMatrix = zeros(n, n+m); % [\lambdaI - A, B]
        testMatrix(:,1:n) = eVal*eye(n) - A;
        testMatrix(:,n+1:end) = B;
        % If rank(testMatrix)<n, then isStable = 0 and system is not stabilizable
        isStable = (rank(testMatrix)==n)*isStable;
    end
end