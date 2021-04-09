function [xk,Pk,yk,Pyy,Pxy] = unscented_transform(X,W_m,W_c,f,h,meas_dim,T)
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
end