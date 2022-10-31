function [xk,Pk,yk,Pyy,Pxy] = unscented_transform(x_k_k,P_k_k,W_m,W_c,f,h,Q,R,meas_dim,T)
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %   author: KURTULUS KEREM SAHİN
    %   e-mail: kurtuluskeremsahin@gmail.com
    %     date: 31.10.2022
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %  This function computes the sufficient statistics of a gaussian
    %  density under a nonlinear transformation using the Unscented
    %  Transform Approach. The state space model used in this model can be
    %  any nonlinear state space model with additive white Gaussian noise.
    %  
    %  State space representation of the model can be given as
    %                x_{k+1} = f(x_k) + w_k 
    %                    y_k = h(x_k) + v_k   
    %   where;
    %            x_k : State vector (in R^n)
    %            y_k : Measurement vector (in R^m)
    %     f:R^n->R^n : State transition function (Possibly Nonlinear)
    %     h:R^n->R^n : Measurement function (Possibly Nonlinear)
    %            w_k : Zero mean white process noise with covariance Q
    %            v_k : Zero mean white measurement noise with covariance R    
    %   INPUT ARGUMENTS
    %          x_k_k : (nx(2n+1) array) Updated state at time k-1
    %          P_k_k : (nxn array) Updated State Covariance Matrix
    %            W_m : (1x(2n+1) array) Weights for mean sigma points 
    %            W_c : (1x(2n+1) array) Weights for covariance sigma points
    %              f : (function handle) Function handle for the state
    %                  transition function
    %              h : (function handle) Function handle for the
    %                  measuremen function
    %              Q : (nxn array) Covariance matrix of the additive white
    %                  process noise
    %              R : (nxn array) Covariance matrix of the additive white
    %                  measurement noise
    %       meas_dim : (scalar) Dimension of the measurement vector
    %              T : (scalar) Sampling Period of the system
    %   OUTPUT ARGUMENTS
    %             xk : (nx1 array) Transformed state mean
    %             Pk : (nxn array) Transformed state covariance matrix
    %             yk : (meas_dimx1 array) Measurement mean vector
    %            Pyy : (meas_dim x meas_dim array) Measurement covariance
    %                  matrix
    %            Pxy : (meas_dim x meas_dim array) State-Measurement cross
    %                  covariance matrix
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    
    % Generation of sigma points
    n = length(x_k_k); 
    X = [x_k_k, x_k_k + kron([1,-1],chol((n/(1-W_m(1)))*P_k_k))];

    % Space allocation and zero initialization for summations
    X_f = zeros(size(X));
    Y_h = zeros(meas_dim, size(X,2));
    xk = zeros(size(X,1),1);
    yk = zeros(meas_dim,1);
    Pk = zeros(size(X,1));
    Pyy = zeros(meas_dim);
    Pxy = zeros(size(X,1),meas_dim);

    % Calculation of the transformed Sigma point Locations and transformed
    % mean values for both state and measurement
    for i = 1:size(X,2)
        X_f(:,i) = f(X(:,i),T);
        Y_h(:,i) = h(X(:,i));
        xk = xk + W_m(i)*X_f(:,i);
        yk = yk + W_m(i)*Y_h(:,i);
    end
    % Covariance Calculations
    for i= 1:size(X,2)
        Pk = Pk + W_c(i)* (X_f(:,i)-xk)*(X_f(:,i)-xk).';
        Pyy = Pyy + W_c(i)* (Y_h(:,i)-yk)*(Y_h(:,i)-yk).';
        Pxy = Pxy + W_c(i)* (X_f(:,i)-xk)*(Y_h(:,i)-yk).';
    end
    Pk = Pk+Q;
    Pyy = Pyy+R;
end