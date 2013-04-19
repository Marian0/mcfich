function q = Q(x)
    q = -1;
endfunction



function [T] = ej1a(dx)
    dt = (0.5*(dx^2))/2;
    disp(dt);
    disp(dx);
    n = 1/dt;
    m = (1/dx)+1;
    T = zeros(n,m);
    for j = 1:m
        T(1,j) = 1 - dx*(j-1); // t = 1 - x
        //T(1,j) = 0;
    end
    for i = 1:n
        T(i,1) = 1;
        T(i,m) = 0;
    end

    for i = 2:n
        for j=2:(m-1)
            T(i,j) = T(i-1,j) + dt/(dx^2)*(T(i-1,j+1) - 2*T(i-1,j) + T(i-1,j-1)) + dt*Q(dx*j);
        end
        plot(T(i,:));    
    end
    

endfunction

    