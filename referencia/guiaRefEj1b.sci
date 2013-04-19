//dx : separacion de nodos
// eps : criterio de parada de error relativo entre iteraciones
//T : campo de temperatura en el dominio
//it : cantidad de iteraciones

//funcion fuente, depende de x
function g = G(x)
    if(x <= 1/2)
        g = 1;
        return;
    end
    g = 0;
endfunction

function [T,it] = ej1b(dx,xl,eps)
    dt = (1/2)*(dx^2) / 2; // Fourier = 0.25
    n = 1/dt;
    m = xl/dx+1;
    T = zeros(n,m);
    x = (0:dx:xl);

    for j = 1:m //condicion inicial
        T(1,j) = 1 - dx*(j-1);
    end
    for i = 1:n //condicion de contorno
        T(i,1) = 1;
    end
    for i = 2:n
        for j = 2:m
            if(j == m) //aplica condicion neumann si estoy en un borde
                value = T(i-1,j-1);
            else 
                value = T(i-1,j+1);
            end
            T(i,j) = dt*G(x(j) - T(i-1,j)) + (dt/(dx^2))*(T(i-1,j-1) - 2*T(i-1,j) + value) + T(i-1,j);
        end
        if( norm( (T(i,:) - T(i-1,:)) ) / norm(T(i,:)) < eps) then
            T = T(1:i,:);
            plot(T($,:));
            it = i;
            return;
        end
    end
    it = floor(n);
    
endfunction

    