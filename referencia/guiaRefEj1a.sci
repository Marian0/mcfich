//dx : separacion de nodos
// eps : criterio de parada de error relativo entre iteraciones
//T : campo de temperatura en el dominio
//it : cantidad de iteraciones


//funcion fuente, testing only
function q = Q(x)
    q = -1;
endfunction

function [T,it] = ej1a(dx,eps)
    dt = (0.5*(dx^2))/2; //se elige un dt para hacer el Fourier = 0.25
    n = 1/dt;
    m = (1/dx)+1;
    T = zeros(n,m);
    for j = 1:m //aplica condiciones iniciales
        T(1,j) = 1 - dx*(j-1); // t = 1 - x
    end
    for i = 1:n //condiciones de contorno
        T(i,1) = 1;
        T(i,m) = 0;
    end
    for i = 2:n 
        for j=2:(m-1)
            T(i,j) = T(i-1,j) + dt/(dx^2)*(T(i-1,j+1) - 2*T(i-1,j) + T(i-1,j-1));
        end
//        plot(T(i,:));    
        if( norm( (T(i,:) - T(i-1,:)) ) / norm(T(i,:)) < eps) then
            T = T(1:i,:);
            it = i;
            return;
        end
        
    end
    it = floor(n);
endfunction

    