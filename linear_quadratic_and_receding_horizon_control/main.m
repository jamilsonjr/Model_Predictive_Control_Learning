clear 
clc
close all
%% Linear Quadratic and Receding Horizon Control Assignment
A = [1.0025 0.1001; 0.0500 1.0025];
B = [0.0050; 0.1001];
C = [1 0];
Q = [5 0; 0 1];
Pf = Q;
R = 0.05;
D = [];
sys = ss(A, B, C, D);

% Discover the shortest value of N that leads to a marginally stable system
for N=1:5
   [K_0, ~] = DP(A, B, N, Q, R, Pf);
   e = eig(A + B*K_0);
   if all( abs(e) < 1 ) 
       break;
   end
end

% Find stationary solution Pf of the Raccati equation. From this Pf, find
% repeat the last step: Whats the difference?

[Pf, ~, ~] = idare(A,B,Q,R);
for N_3=1:5
   [K_0, ~] = DP(A, B, N_3, Q, R, Pf);
   e = eig(A + B*K_0);
   if all( abs(e) < 1 ) 
       break;
   end
end