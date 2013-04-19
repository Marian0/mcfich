//SALIDA:
// x : vector solucion Phi
// A : Matriz del problema (stencil)
// b : vector de terminos independientes del stencil
// phi : solucion Phi en forma matricial
// u : campo de velocidades en x
// v : campo de velocidades en y

function [x,A,b,phi,u,v] = ej4(x0,xl,y0,yl,dx) //supone dx = dy
    n = (((xl-x0)/dx)+1)*(((yl-y0)/dx)+1);
    A = zeros(n,n);
    alfa = (xl-x0)/dx+1;
    b = zeros(n,1);
    xvec = [x0:dx:xl];
    for k = 1:n
        A(k,k) = -4;
        //pruebo las esquinas primero
        if(k == 1) then //abajo izquierda
            A(k,k+1) = 2;
            A(k,k+alfa) = 2;
        elseif(k == alfa) then //arriba izquierda
            A(k,k+alfa) = 2;
            A(k,k-1) = 2;
            b(k) = -xvec(k/alfa);    
        elseif(k == (n-alfa+1)) then //abajo derecha
            A(k,k+1) = 2;
            A(k,k-alfa) = 2;
        elseif(k == n) then //arriba derecha
            A(k,k-alfa) = 2;
            A(k,k-1) = 2;
            b(k) = -xvec(k/alfa);
        elseif((k < alfa) & (k <> 1)) then    //pruebo los bordes, borde izquierdo
                A(k,k+alfa) = 2;
                A(k,k+1) = 1;
                A(k,k-1) = 1;
        elseif((k > (n-alfa+1)) & (k < n)) then //borde derecho
                    A(k,k-alfa) = 2;
                    A(k,k+1) = 1;
                    A(k,k-1) = 1;
        elseif(modulo(k,alfa) == 1) then //borde inferior
                A(k,k-alfa) = 1;
                A(k,k+alfa) = 1;
                A(k,k+1) = 2;
        elseif((modulo(k,alfa) == 0) & (k <> alfa) & (k <> n)) then //borde superior
                A(k,k-alfa) = 1;
                A(k,k+alfa) = 1;
                A(k,k-1) = 2;
                b(k) = -xvec(k/alfa);
        else //nodos interiores
            A(k,k+1) = 1;
            A(k,k-1) = 1;
            A(k,k-alfa) = 1;
            A(k,k+alfa) = 1;
        end
    end
    x = A\b;
    m = n/alfa;
    phi = matrix(x,[m,m]); //convierte a matriz el vector solucion
    u = zeros(m,m); //campo de velocidades en x
    v = zeros(m,m); //campo de velocidades en y
    u(1:m,1) = 0; //primer condicion Neumann
    u(1:m,m) = 0; //primer condicion Neumann
    v(m,1:m) = 0; //segunda condicion Neumann
    v(1,1:m) = xvec(1:m); //tercera condicion Neumann
    for i=2:(m-1)
        for j=2:(m-1)
            u(i,j) = (phi(i+1,j) - phi(i-1,j))/(2*dx);
            v(i,j) = (phi(i,j+1) - phi(i,j-1))/(2*dx);
        end
    end
endfunction
