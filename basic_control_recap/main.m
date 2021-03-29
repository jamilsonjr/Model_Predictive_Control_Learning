clear
close all
%% First assignment
%% Question 1
% Matrices
Ac = [ 0    1;
      0.95 0];
Bc = [0  ;
     1 ];    
Cc = [ 1 0 ]; 
Dc = [];
h = 0.1;
%% System
sys_cont = ss(Ac, Bc, Cc, Dc);
sys_disc = c2d(sys_cont, h);
Ad = sys_disc.a;
Bd = sys_disc.b;
Cd = sys_disc.c;
Dd = sys_disc.d;
%% Question 2
tau = 0.5 * h;
syms s;
A_ = expm(Ac*h);
B_1 = double( expm( Ac*(h-tau) )* int( expm( Ac*s )*Bc,s,0, tau ) );
B_2 = double( int( expm( Ac*s )*Bc, s, 0, h-tau ) );

Aa = [A_ B_1; zeros(1,2) 0];
Ba = [B_2; 1];
Ca = [Cc 0];
Da = zeros();
sys_aug = ss(Aa, Ba, Ca, Da); 

%% Question 3
% Controlability
Co_cont = rank( ctrb( sys_cont ) );
Co_aug = rank( ctrb( sys_aug ) );
Co_disc = rank( ctrb( sys_disc ) );

% Observability
Obs_cont = rank( obsv( sys_cont ) );
Obs_aug = rank( obsv(sys_aug) );
Obs_disc = rank( obsv(sys_disc) );

%% Question 6
lambda_1 = -4 + 6*1j;
lambda_2 = -4 - 6*1j;

pole_1 = exp( h* ( -4 + 6*1j ) );
pole_2 = exp( h* ( -4 - 6*1j ) );

cl_poles = [pole_1; pole_2];

K_disc = place(Ad, Bd, cl_poles);
K_aug1 = [K_disc 0];
K_aug2 = place(Aa, Ba, [cl_poles; 0]);

%have labels and legend.
cl_sys_d = ss( Ad-Bd*K_disc, Bd, Cd, Dd, h );
cl_sys_a1 = ss( Aa-Ba*K_aug1, Ba, Ca, Da, h );% Closed loop system using feedback from system (3).
cl_sys_a2 = ss( Aa-Ba*K_aug2, Ba , Ca, Da, h );% Closed loop system with new pole placement.
figure();
hold on;
step(cl_sys_d)
step(cl_sys_a1)
legend('System (3)','System (4)');
figure();
hold on;
step(cl_sys_d)
step(cl_sys_a1)
step(cl_sys_a2)
legend('System (3)','System (4)','System (4) with New Feedback');
%% Question 7
ys = pi/6;
steady_state_matrix = [(eye(length(Aa))-Aa) -Ba; Ca 0];
x_steady_state = steady_state_matrix \ [zeros(3,1); ys];

% reference gain
K_reference = ( Ca * ( (eye(3) - Aa + Ba*K_aug2 )^(-1) )*Ba )^(-1);

N=50;
t=0:h:(N-1)*h;
x=zeros(3,N);
y=zeros(1,N);
x(:,1)=[0;0;0];
y(:,1)=0;
for k=1:N-1
    x(:,k+1) = (Aa-Ba*K_aug2)*x(:,k) + Ba*K_reference * ys;
    y(:,k+1) = Ca*x(:,k+1);
end
figure();  
plot(t,y);
xlabel('Time(Seconds)');
ylabel('Amplitude');
title('Step Response of System (4) with Reference Tracking');
grid on

%% Question 8
Bd = [0 1 0]';
Ae = [Aa Bd; zeros(1,length(Aa)) eye(1)];
Be = [Ba' 0]';
Ce = [Ca 0];
De = zeros();

sys_dist = ss(Ae, Be, Ce, De, h);

% Controlability
Co_e = rank( ctrb( sys_dist ) ); % Not controllable
Obs_e = rank( obsv( sys_dist ) ); % Observable

poles_unreachable_state = eig(Ae);
PBH = zeros(1,length(poles_unreachable_state));
for i=1:length(poles_unreachable_state)
    PBH(i) = rank( [poles_unreachable_state(i)*eye(length(Ae))-Ae Be] );
end
