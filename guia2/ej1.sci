function fi = simp_compfk(a,b,n,l,m)
    h = (b-a)/n;
    x = zeros(n+1,1);
    for i=0:n
        x(i+1) = a+i*h;
    end

    suma1 = 0;
    for j=1:((n/2)-1)
        suma1 = suma1 + fk(x(2*j),l,m);
    end
    
    suma2 = 0;
    for j=2:(n/2)
        suma2 = suma2 + fk(x(2*j-1)),l,m;
    end
    
    fi = h/3*(fk(x(1),l,m) + 2*suma1 + 4*suma2 + fk(x(n),l,m));
endfunction

    
    


function fi = simp_compff(a,b,n,l)
    h = (b-a)/n;
    x = zeros(n+1,1);
    for i=0:n
        x(i+1) = a+i*h;
    end

    suma1 = 0;
    for j=1:((n/2)-1)
        suma1 = suma1 + F(x(2*j),l);
    end
    
    suma2 = 0;
    for j=2:(n/2)
        suma2 = suma2 + F(x(2*j-1),l);
    end
    
    fi = h/3*(F(x(1),l) + 2*suma1 + 4*suma2 + F(x(n),l));
endfunction
        
    
    


function z = simpF(a,b,l)
    la = (a+b)/2;
    z = (b-a)/6*(F(a,l) + F(la,l) + F(b,l));
endfunction

function z =simpK(a,b,l,m)
    z = (b-a)/6*(fk(a,l,m) + fk((a+b)/2,l,m) + fk(b,l,m));
endfunction

function z = F(x,l)
    z = x^l*(1-x)*(sin(%pi*x/2) - x);
endfunction

function z = fk(x,l,m)
    z = (x^l)*(x^m)*(1-x)^2;
endfunction
    

function [a,K,f] = guia2ej1(a,b,dx)
    l = (b-a)/dx ;
    a = zeros(l,1);
    f = a;
    K = zeros(l,l);
    for i=1:l
        for j=1:l
            K(i,j) = simp_compfk(0,1,10,i,j);
        end
        f(i) = simp_compff(0,1,10,i);
    end
    a = K\f;
endfunction


function y = final(x,x0)
    n = length(x);
    y = 0;
    for i=1:n
        y = y+(1+x0)+x(i)*(x0^i * (1-x0));
    end
endfunction


function [t,v] = plotfinal(x,a,b,dx)
    n = (b-a)/dx;
    t = zeros(n,1);
    v = zeros(n,1);
    for i=1:n
        t(i) = final(x,(i-1)*dx);
        v(i) = phi((i-1)*dx);
    end
endfunction

function y = phi(x)
    y = 1+sin(%pi*x/2);
endfunction
