% clear all
% close all 
% path_genreration
sample_count = 400;

meas_dim = 2; 
T = 1 ;
n = 5;
sigma_state_filter = 20;
sigma_meas_filter = 40;

W_0 = -0.3;

Q = [sigma_state_filter^2  0
            0              0.1];
        
R = sigma_meas_filter^2  * eye(2);


f = @(mu,t)state_transition_cartesian(mu,t);
H = @(mu) h(mu);
x = zeros(n,size(path_org,2));
x_0 = [0 ; 0; 10 ;10 ; 1e-3 ];
x(:,1) = x_0;
P_k = 10 * eye(n);

for k = 1: sample_count-1
    %% Find Sigma Points X
    X = zeros(n,2*n+1);
    X(:,1) = x(:,k) ;
    Sigma_sqrt = chol((n/(1-W_0))*P_k);
    X(:,2:n+1) = X(:,1) - Sigma_sqrt ;
    X(:,n+2:end) = X(:,1) + Sigma_sqrt ;
    W_m = [W_0; ((1-W_0)/(2*n).*ones(2*n,1))];
    W_c = W_m ;
    
    %% Time Update
    [xk_hat,Pk_hat,yk_hat,Pyy,Pxy] = unscented_transform(X,W_m,W_c,f,H,meas_dim,T);
    theta = atan2(xk_hat(4),xk_hat(3));
    B = [(0.5*T^2*cos(theta))       0
         (0.5*T^2*sin(theta))       0
               T*cos(theta)         0
               T*sin(theta)         0
                    0               T    ];
    cov_w = B* Q* B.';
    
    Pk_hat = Pk_hat + cov_w;
    Pyy = Pyy + R;

    %% Measurement Update
    Kk = Pxy / Pyy ;
    x(:,k+1) = xk_hat + Kk * (measurement(:,k) - yk_hat);
    P_k = Pk_hat - Kk * Pyy * Kk.';
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




function z = h(x)
    z= [1 0 0 0 0;0 1 0 0 0 ]*x;
end