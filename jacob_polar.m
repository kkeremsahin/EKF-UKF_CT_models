function [x_k1,F] =  jacob_polar(x,T)
x1 = x(1);
x2 = x(2);
v = x(3);
h = x(4);
w = x(5);

x_k1 = [x1 + (2*v*cos(h + (T*w)/2)*sin((T*w)/2))/w
        x2 + (2*v*sin(h + (T*w)/2)*sin((T*w)/2))/w
                              v
                           h + T*w
                              w                     ];
                          
F = [ 1, 0, (2*cos(h + (T*w)/2)*sin((T*w)/2))/w, -(2*v*sin(h + (T*w)/2)*sin((T*w)/2))/w, (T*v*cos(h + (T*w)/2)*cos((T*w)/2))/w - (T*v*sin(h + (T*w)/2)*sin((T*w)/2))/w - (2*v*cos(h + (T*w)/2)*sin((T*w)/2))/w^2
     0, 1, (2*sin(h + (T*w)/2)*sin((T*w)/2))/w,  (2*v*cos(h + (T*w)/2)*sin((T*w)/2))/w, (T*v*cos(h + (T*w)/2)*sin((T*w)/2))/w - (2*v*sin(h + (T*w)/2)*sin((T*w)/2))/w^2 + (T*v*sin(h + (T*w)/2)*cos((T*w)/2))/w
    0, 0,                                   1,                                      0,                                                                                                                       0
    0, 0,                                   0,                                      1,                                                                                                                       T
    0, 0,                                   0,                                      0,                                                                                                                       1];
end
                          