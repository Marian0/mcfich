function [x] = armarQpotencial(lx,ly,dx,dy)
    n = (lx/dx+1)*(ly/dy+1);
    x = zeros(n,1);
    x(1) = 10;
    //x(n) = 0;
endfunction


function [x,A,b] = potencial(lx,ly,dx,dy)
    I = lx/dx+1;
    J = ly/dy+1;
    n = (I)*(J);
    A = zeros(n,n);
    x = zeros(n,1);
    for i=1:n
        A(i,i) = -2*(1/(dx^2)+1/(dy^2));
        if(i > 1) then
            A(i,i-1) = 1/(dy^2);
        end
        if(i < n) then
            A(i,i+1) = 1/(dy^2);
        end
        if((i-4) >= 1) then
            A(i,i-4) = 1/(dx^2);
        end
        if((i+4) <= n) then
            A(i,i+4) = 1/(dx^2);    
        end
    end
    b = armarQpotencial(lx,ly,dx,dy);
    x = A\b;
endfunction
