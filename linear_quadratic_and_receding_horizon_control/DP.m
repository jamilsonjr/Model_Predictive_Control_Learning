function [K_0,P_0]=DP(A,B,N,Q,R,Pf)
    % Obtain the last values for K and P
    K = -((R + B'*Pf*B)^(-1)) * B'*Pf*A;
    P = Q + A'*Pf*A - A'*Pf*B*(( R + B'*Pf*B )^(-1)) * B'*Pf*A;
    % Obtain the all the way down to the first value
    for i=N-1:-1:1
        K = -((R + B'*P*B)^(-1)) * B'*P*A;
        P = Q + A'*P*A - A'*P*B*(( R + B'*P*B )^(-1)) * B'*P*A;
    end
    K_0 = K;
    P_0 = P;
end

