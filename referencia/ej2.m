
%Calcula la matriz B para elementos impares
%  Bodd
%-1    0  -1   0  0   0
% 0  -dx   0   0  0   dx
%-dx  -1   0  -1  dx  0
function B = Bodd(dx)
    B = zeros(3,6);
    B = [-1 0 1 0 0 0; 0 -2 0 0 0 2; -2 -1 0 1 2 0];
    area = dx/2;
    B = 1/(2*area) * B;
endfunction

%Calcula la matriz B para elementos pares
%  Beven
%  0    0   1   0  -1   0
%  0  -dx   0  dx   0   0
%-dx    0  dx   1   0  -1
function B = Beven(dx)
    B = zeros(3,6);
    B = [0 0 1 0 -1 0; 0 -dx 0 dx 0 0;-dx 0 dx 1 0 -1];
    area = dx/2;
    B = 1/(2*area) * B;
endfunction

%Calcula la matriz de conectividades:
%C
%1 2 3
%2 3 4
%...
%n-2 n-1 n
function C = conex(cantElementos)
    C = zeros(cantElementos,3);
    for i=1:cantElementos
        C(i,1) = i;
        C(i,2) = i+1;
        C(i,3) = i+2;
    end
endfunction

%Calcula la matriz D
%D
%E/(1-nu*nu) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2]

function D = Dmatrix(E,nu)
    D = [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    D = (E/(1-nu*nu)) * D;
endfunction

%Calcula la matriz K
%n: numero de elementos (triangulos)
%   Todos van a tener altura 1, y base = 10/(n/2)
%t: espesor
%E: Modulo de Young
%nu: Coeficiente de Poisson
%  Forma de la malla:
% 2___4
% |\  |
% | \ |  <-1
% |__\|
% 1   3
%   ^ 
%   dx
%
% Elemento 1: [1 2 3]
% Elemento 2: [2 3 4]
function K = Kmatrix(n,t,E,nu)
    %dx es el lado mas ancho del triangulo
    %dy siempre vale 1
    dx = 10/(n/2);
    %Matrices B para los 2 tipos de elementos.
    Bimpar = Bodd(dx);
    Bpar = Beven(dx);
    %Matriz de conectividades
    C = conex(n);
    D = Dmatrix(E,nu);
    %Calcula las K elementales
    areaElemento = dx/2;
    KeImpar = Bimpar' * D * Bimpar * t * areaElemento;
    KePar = Bpar' * D * Bpar * t * areaElemento;
    
    K = zeros((n+2)*2,(n+2)*2);
    for i=1:n
        fila1 = 2*C(i,1)-1;
        fila2 = 2*C(i,2)-1;
        fila3 = 2*C(i,3)-1;
        if(mod(i,2) == 1) %impar
           K([fila1,fila1+1],[fila1,fila1+1]) = K([fila1,fila1+1],[fila1,fila1+1]) + KeImpar(1:2,1:2);
           K([fila1,fila1+1],[fila2,fila2+1]) = K([fila1,fila1+1],[fila2,fila2+1]) + KeImpar(1:2,3:4);
           K([fila1,fila1+1],[fila3,fila3+1]) = K([fila1,fila1+1],[fila3,fila3+1]) + KeImpar(1:2,5:6);
           K([fila2,fila2+1],[fila1,fila1+1]) = K([fila2,fila2+1],[fila1,fila1+1]) + KeImpar(3:4,1:2);
           K([fila2,fila2+1],[fila2,fila2+1]) = K([fila2,fila2+1],[fila2,fila2+1]) + KeImpar(3:4,3:4);
           K([fila2,fila2+1],[fila3,fila3+1]) = K([fila2,fila2+1],[fila3,fila3+1]) + KeImpar(3:4,5:6);
           K([fila3,fila3+1],[fila1,fila1+1]) = K([fila3,fila3+1],[fila1,fila1+1]) + KeImpar(5:6,1:2);
           K([fila3,fila3+1],[fila2,fila2+1]) = K([fila3,fila3+1],[fila2,fila2+1]) + KeImpar(5:6,3:4);
           K([fila3,fila3+1],[fila3,fila3+1]) = K([fila3,fila3+1],[fila3,fila3+1]) + KeImpar(5:6,5:6);
       else %par
           K([fila1,fila1+1],[fila1,fila1+1]) = K([fila1,fila1+1],[fila1,fila1+1]) + KePar(1:2,1:2);
           K([fila1,fila1+1],[fila2,fila2+1]) = K([fila1,fila1+1],[fila2,fila2+1]) + KePar(1:2,3:4);
           K([fila1,fila1+1],[fila3,fila3+1]) = K([fila1,fila1+1],[fila3,fila3+1]) + KePar(1:2,5:6);
           K([fila2,fila2+1],[fila1,fila1+1]) = K([fila2,fila2+1],[fila1,fila1+1]) + KePar(3:4,1:2);
           K([fila2,fila2+1],[fila2,fila2+1]) = K([fila2,fila2+1],[fila2,fila2+1]) + KePar(3:4,3:4);
           K([fila2,fila2+1],[fila3,fila3+1]) = K([fila2,fila2+1],[fila3,fila3+1]) + KePar(3:4,5:6);
           K([fila3,fila3+1],[fila1,fila1+1]) = K([fila3,fila3+1],[fila1,fila1+1]) + KePar(5:6,1:2);
           K([fila3,fila3+1],[fila2,fila2+1]) = K([fila3,fila3+1],[fila2,fila2+1]) + KePar(5:6,3:4);
           K([fila3,fila3+1],[fila3,fila3+1]) = K([fila3,fila3+1],[fila3,fila3+1]) + KePar(5:6,5:6);
        end    
    end
endfunction

%Calcula la fuerza elemental.
% Los elementos impares corresponden a la direccion x
% Los elementos pares son los correspondientes a la direccion y
function F = Felem(px,py)
    F = [0 0.5 0 0.5 0 0]';
    F(1:2:6) = F(1:2:6)*px;
    F(2:2:6) = F(2:2:6)*py;
endfunction

%Calcula el vector F completo
%La fuerza solo se aplica en la cara vertical del ultimo elemento
% Entonces, el vector F es cero EXCEPTO en los ultimos 6 valores, que aplica el vector Felem
function F = Fvector(n,px,py)
    F = zeros((n+2)*2,1);
    Fe = Felem(px,py);
    C = conex(n);
    idx = (n+2)*2;
    F(idx-5:idx,1) = Fe;
endfunction

function [UV,K,F] = resolver(n,t,E,nu,px,py)
    F = Fvector(n,px,py);
    K = Kmatrix(n,t,E,nu);
    ultimoidx = 2*(n+2);
    %Aplica condicion que los nodos 1 y 2 estan fijos, borro filas de 1 a 4
    K = K(5:ultimoidx,5:ultimoidx);
    F = F(5:ultimoidx,1);
    %resuelve
    UV = K\F;
    UV = [0;0;0;0;UV]; %agrego los dos nodos del borde a la solucion
endfunction
    
    

    