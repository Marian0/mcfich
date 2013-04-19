function y = qfunction(x)
    if(x <= 0.5) then
        y = 1;
        return;
    end
    y = 0;
endfunction

    

function [A,b,x] = matriz(h,c0,cN)
    A = zeros(1/h - 1,1/h - 1);
    n = size(A,'r');
    q = zeros(n,1);
    for i=1:n
        q(i) = qfunction(i*h);
        
    end 
    
    //disp(q)   
    b = q;
    b = b*-(h^2);
    
    b(1) = b(1) - c0;
    b(n) = b(n) - cN;
    
    for i=1:n
        for j=1:n
            if(i == j) then
                A(i,j) = -2;
            end
            if((i == j-1) | (i == j+1)) then
                A(i,j) = 1;
            end
        end
    end
    
    x = A\b;
endfunction
