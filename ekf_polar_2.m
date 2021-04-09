% clearvars; close all;
% This is the coordinated turn model extended kalman filter for 2D tracking problem
% with states of the form x = [x1,x2,v,h,w]^T and the measurements of the
% form y = [x1,x2]^T


T = 1; % Sampling Time (s)
sigma_state_filter = 50;
sigma_meas_filter = 20;
%% State Matrix Definitions
% A is the jacobian of the state transition function f(x(k))
% B is chosen from two alternatives as G_p1
B = [0       0
     0       0
     T       0
     0       0
     0       T  ];
    
C = [ 1 0 0 0 0
      0 1 0 0 0];
Q = [sigma_state_filter^2     0
            0                 0.1];
R = sigma_meas_filter^2  * eye(2);

%% Path and measurement generation
% Path Parameters
sample_count = 400;
v= 10;
sigma_a = 0.2;
sigma_alpha = pi/200;
w_1 = pi/200;
w_2 = pi/200;
% Path Generation
%[path_org, vel_org, w_org] = generate_turn_path(T,sample_count,v,w_1,w_2,sigma_a,sigma_alpha);
% Measurement generation;
sigma_meas = 20 ;
% measurement = path_org + sigma_meas * randn(size(path_org));

%% Linear Kalman Filter
% Initial State definitions and space allocation for states
x = zeros(5,size(path_org,2));
x_0 = [0 ; 0; 10; 1; 1e-3];
x(:,1) = x_0;
P_k = 20 * eye(5);

for k = 1: sample_count-1
    x_k = x(:,k);
    
    % Time_update
    
    [x_kp1_hat, A] = jacob_polar(x_k,T);
    B = [0.5*T^2*cos(x_kp1_hat(4))       0
         0.5*T^2*sin(x_kp1_hat(4))       0
         T                            0
         0                          0.5*T^2
         0                            T  ];
    
    P_k_hat = A*P_k* A.' + B* Q* B.';
    
    % Measurement update
    S=C*P_k_hat*C'+R;
    
    Pk=P_k_hat-P_k_hat*C'/S*C*P_k_hat;
    
    x(:,k+1)=x_kp1_hat + Pk*C'/S*(measurement(:,k)-C*x_kp1_hat);
    
end

%% Plot the Results
figure
plot (path_org(1,:),path_org(2,:),'Linewidth',2)
hold on
plot(x(1,:),x(2,:),'Linewidth',2)
scatter(measurement(1,:),measurement(2,:))
%% Find Errors
Error = path_org - [x(1,:);x(2,:)];
rmse = sqrt(Error(1,:).^2 + Error(2,:).^2);
figure
plot(rmse);




    