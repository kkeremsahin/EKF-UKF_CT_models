clearvars; close all;

load ('ex_ground_truth.mat')
seed = 1000;
rng(seed);
T = 1; % Sampling Time (s) 
ground_truth_state = [position_gt;velocity_gt;omega_gt];
ground_truth_state = ground_truth_state(:,1:(fs*T):end);  %downsample the generated path
N_samples = size(ground_truth_state,2);
sigma_state_filter = 10;
sigma_meas_filter = 10;

meas_dim = 2; 
n = 5;

W_0 = -0.3;

Q = [sigma_state_filter^2  0
           0            0.1];
C = [eye(2),zeros(2,3)];
R = sigma_meas_filter^2  * eye(2);
B = [0       0
     0       0
     T       0
     0       0
     0       T  ];

cov_w = B* Q* B.';
f = @(mu,t)state_func_polar(mu,t);
H = @(mu) h(C,mu);

%% Measurement generation
measurement = C* ground_truth_state + mvnrnd(zeros(size(C,1),1),R,N_samples).';

%% Unscented Filtering
x = zeros(n,size(ground_truth_state,2));
x_0 = [0 ; 0; 27; 0; 1e-4];
x(:,1) = x_0;
P_k = 10 * eye(n);
W_m = [W_0; ((1-W_0)/(2*n).*ones(2*n,1))];
W_c = W_m ;
for k = 1: N_samples-1
    %% Time Update
    [xk_hat,Pk_hat,yk_hat,Pyy,Pxy] = unscented_transform(x(:,k),P_k,W_m,W_c,f,H,cov_w,R,meas_dim,T);

    %% Measurement Update
    Kk = Pxy / Pyy ;
    x(:,k+1) = xk_hat + Kk * (measurement(:,k) - yk_hat);
    P_k = Pk_hat - Kk * Pyy * Kk.';
end

%% Plot the Results
figure
plot (ground_truth_state(1,:),ground_truth_state(2,:),'Linewidth',2)
hold on
plot(x(1,:),x(2,:),'Linewidth',2)
scatter(measurement(1,:),measurement(2,:))

function z = h(C,x)
    z= C*x;
end