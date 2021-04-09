% clearvars; close all;
% This is the coordinated turn model extended kalman filter for 2D tracking problem
% with states of the form x = [x1,x2,x1',x2',w]^T and the measurements of the
% form y = [x1,x2]^T

% The constant velocity model is implemented;
% x(k+1) = f(x(k)) + G(x(k))(v(k))
% y(k) = C x(k) + w(k)
% v(k)~ N(0,Q)
% w(k)~ N(0,R)
% f(x(k)) = [ x1 + (v1/w)*sin(w*T) - (v2/w)*(1-cos(w*T))
%             x2 + (v1/w)*(1-cos(w*T)) - (v2/w)*sin(w*T)
%                      v1*cos(w*T) - v2*sin(w*T)
%                      v1*sin(w*T) + v2*cos(w*T)      
%                                    w                  ] 

T = 1; % Sampling Time (s)
sigma_state_filter = 10;
sigma_meas_filter = 10;
%% State Matrix Definitions
% A is the jacobian of the state transition function f(x(k))
% B is chosen from two alternatives as G_c2

    
C = [ 1 0 0 0 0
      0 1 0 0 0];
Q = [sigma_state_filter^2  0
            0             0.1];
R = sigma_meas_filter^2  * eye(2);


sample_count = 400;


%% Linear Kalman Filter
% Initial State definitions and space allocation for states
x = zeros(5,size(path_org,2));
x_0 = [0 ; 0; 10; 10; 1e-3];
x(:,1) = x_0;
P_k = 20 * eye(5);

for k = 1: sample_count-1
    x_k = x(:,k);
    theta = atan2(x_k(4),x_k(3));
    
    B = [(0.5*T^2*cos(theta))       0
         (0.5*T^2*sin(theta))       0
               T*cos(theta)         0
               T*sin(theta)         0
                    0               T    ];
 
    % Time_update
    
    [x_kp1_hat, A] = jacob_cartesian(x_k,T);
    
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




    