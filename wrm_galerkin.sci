function z = x(l,dx)
    z = (l-1)*dx;
endfunction

function z = y(m,dy)
    z = (m-1)*dy;
endfunction

function z = func(l,m,x,y)
    z = wlm(l,m,x,y) * nlm(l,m,x,y);
endfunction

//se define wlm = nlm por ser Galerkin
function z = wlm(l,m,x,y)
    z = nlm(l,m,x,y);
endfunction

function z = nlm(l,m,x,y)
    z = sin(l*%pi*torad(x))*sin(m*%pi*torad(y));
endfunction

function rad = torad(g)
    rad = (%pi/180) * g;
endfunction


function z = trap(l,m,dx,dy)
    a = x(l,dx); //x_l
    b = x(l+1,dx); // x_l+1
    c = y(m,dy);//y_m
    d = y(m+1,dy);//y_m+1
    z = ((d-c)*(b-a)/4) * (func(l,m,a,c) + func(l+1,m,b,c) + func(l,m+1,a,d) + func(l+1,m+1,b,d));
endfunction

function [a,K,f] = wrm_galerkin(x0,xl,y0,yl,dx,dy)
    l = (xl-x0)/dx+1;
    m = (yl-y0)/dy+1;
    K = zeros(l,m);
    a = zeros(l*m,1);
    a = [0,0,0,0,0,0,0.5,0.75,0.25,0,0,1,1.75,0.5,0,0,0.75,1.25,0.25,0,0,0,0,0,0];
    f = a;
    for i = 1:l
        for j = 1:m
            K(i,j) = trap(i,j,dx,dy);
        end
    end
endfunction
