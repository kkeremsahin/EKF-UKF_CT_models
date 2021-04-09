function [x_k1] =  state_func_polar(x,T)
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
end