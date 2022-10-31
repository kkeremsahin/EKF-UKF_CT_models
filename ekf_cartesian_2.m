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

%% State Matrix Definitions
% A is the jacobian of the state transition function f(x(k))
% B is chosen from two alternatives as G_c2

    
C = [eye(2),zeros(2,3)];
Q = [sigma_state_filter^2  0
            0             0.01];
R = sigma_meas_filter^2  * eye(2);


%% Measurement generation
measurement = C* ground_truth_state + mvnrnd(zeros(size(C,1),1),R,N_samples).';

%% Extended Kalman Filter
% Initial State definitions and space allocation for states
x = zeros(5,N_samples-1);
x_kp1_hat = [0 ; 0; 27; 0; 1e-4];
x(:,1) = x_kp1_hat;
P_k_hat = diag([100*ones(1,4),1e-5]);
for k = 1:N_samples-1 
    % Measurement update
    S = C*P_k_hat*C'+R;
    
    Kk = P_k_hat* C.'/S;

    P_k = P_k_hat- Kk*S*Kk.';
    
    x(:,k)= x_kp1_hat + Kk*(measurement(:,k)-C*x_kp1_hat);
  
    % Time_update
    theta = atan2(x(4,k),x(3,k));
    B = [(0.5*T^2*cos(theta))       0
         (0.5*T^2*sin(theta))       0
               T*cos(theta)         0
               T*sin(theta)         0
                    0               T    ];
 
    
    [x_kp1_hat, A] = jacob_cartesian(x(:,k),T);
    
    P_k_hat = A*P_k* A.' + B* Q* B.';
end

%% Plot the Results
figure
plot (ground_truth_state(1,:),ground_truth_state(2,:),'Linewidth',2)
hold on
plot(x(1,:),x(2,:),'Linewidth',2)
plot(measurement(1,:),measurement(2,:),'o')