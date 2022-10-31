function [x_k1,F] =  jacob_cartesian(x,T)
x1 = x(1);
x2 = x(2);
v1 = x(3);
v2 = x(4);
w  = x(5);

x_k1 = [ x1 + (v1/w)*sin(w*T) - (v2/w)*(1-cos(w*T))
        x2 + (v1/w)*(1-cos(w*T)) + (v2/w)*sin(w*T)
        v1*cos(w*T) - v2*sin(w*T)
        v1*sin(w*T) + v2*cos(w*T)
        w                     ] ;

F = [   1, 0,        sin(T*w)/w, (cos(T*w) - 1)/w, (T*v1*cos(T*w))/w - (v2*(cos(T*w) - 1))/w^2 - (v1*sin(T*w))/w^2 - (T*v2*sin(T*w))/w
        0, 1, -(cos(T*w) - 1)/w,       sin(T*w)/w, (v1*(cos(T*w) - 1))/w^2 - (v2*sin(T*w))/w^2 + (T*v2*cos(T*w))/w + (T*v1*sin(T*w))/w
        0, 0,          cos(T*w),        -sin(T*w),                                                     - T*v2*cos(T*w) - T*v1*sin(T*w)
        0, 0,          sin(T*w),         cos(T*w),                                                       T*v1*cos(T*w) - T*v2*sin(T*w)
        0, 0,                 0,                0,                                                                                   1];
end
