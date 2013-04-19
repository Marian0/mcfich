function [Q] = armarQ(lx,ly,dx,dy)
    Q = zeros((lx/dx+1)*(ly/dy+1),1);
    I = lx/dx+1;
    J = ly/dy+1;
    paso = 1;
    for i = 1:(I-1)
        for j = 1:(J-1)
            Q(paso) = -2*((i*dx)^2 + (j*dy)^2) * dx^2 * dy^2;
            paso = paso + 1;
        end
    end
endfunction


function [x,A,Q] = matriz2d(lx,ly,dx,dy)
    n = (lx/dx+1)*(ly/dy+1); //cantidad de incognitas
    I = lx/dx+1;
    J = ly/dy+1;
    A = zeros(n,n);
    Q = armarQ(lx,ly,dx,dy);
   
    for i = 1:(n)
       A(i,i) = -2*(dx^2 + dy^2);
       if((i-1) > 0) then
           A(i,i-1) = dx^2;
       end
       if((i+1) <= n) then
           A(i,i+1) = dx^2;
       end
       if((i-5) > 0) then
           A(i,i-5) = dy^2;
       end
       if((i+5) <= n) then
           A(i,i+5) = dy^2;
       end
    end
    x = A\Q;
endfunction  