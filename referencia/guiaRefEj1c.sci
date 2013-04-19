//A0 : Condicion Dirichlet en x=0
//dx : separacion de nodos
//l0,lx : dimension de la barra


//funcion fuente, depende de x
function f = F(x)
    if(x <= 1)
        f = -10;
    else if(x <= 2)
            f = 5;
        else 
            f = -1;
        end
    end
endfunction


function [x,A,b] = ej1c(A0,dx,l0,lx)
    n = (lx-l0)/dx;
    A = zeros(n,n); //stencil
    A(1,1) = -2;
    A(1,2) = 1;
    for i = 2:(n-1)
        A(i,i) = -2;
        A(i,i-1) = 1;
        A(i,i+1) = 1;
    end
    A(n,n-1) = 2;
    A(n,n) = -2;
    
    b = zeros(n,1);

    pos = [l0:dx:lx];
    for i = 1:n
        b(i) = (dx^2) * (-1) * F(pos(i)); //terminos independientes
    end
    b(1) = b(1) - A0; //por cond. dirichlet
    x = A\b;
    x = [A0;x]; //agrega el valor inicial de Dirichlet
endfunction