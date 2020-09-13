%==============================================================================
% MOD4 4MC00 / Jelle Langedijk / TU/e
% TEMPLATE V1.1 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================

%% PREREQUISITIES
clear all; close all; clear vars; clc;


 

%% VERTICAL LOAD
N = 8;
h = 1/(N);                               
xi = [0:h:N];


%% MATRIX FILLER (Neumann)
A = zeros(N);
A(1,1) = (2)/(3*h^(2));
A(1,2) = (-2)/(3*h^(2));


A(N,N-1) = (-2)/(3*h^(2));
A(N,N) = (2)/(3*h^(2));

for i = 2:N-1
    A(i,i-1) = -(1)/(h^(2));
end   

for j = 2:N-1
    A(j,j) = (2)/(h^(2));
end  

for k = 2:N-1
    A(k,k+1) = -(1)/(h^(2));
end

%% MANUFACTURED SOLUTIONS

%% SINGULARITY

c = cond(A);

%% MATRIX FILLER (Dirichlet)
A = zeros(N);
A(1,1) = (2)/(h^(2));
A(1,2) = (-1)/(h^(2));


A(N,N-1) = (-1)/(h^(2));
A(N,N) = (2)/(h^(2));

for i = 2:N-1
    A(i,i-1) = -(1)/(h^(2));
end   

for j = 2:N-1
    A(j,j) = (2)/(h^(2));
end  

for k = 2:N-1
    A(k,k+1) = -(1)/(h^(2));
end

cd = cond(A);

%% SOLUTION
f = 1;
for i = 1:N
    b(i) = f;
end

b = transpose(b);

U = A\b;

%% LU DECOMPOSITION

%% GAUSS-SEIDEL RELAXATION
stop = 20; 

A_tilde = A;

for i = 1:N
    A_tilde(i,i+1) = 0;
end

A_tilde = A_tilde(:,1:N);

A_tilde_inv = inv(A_tilde);
 
x_k = zeros(N,stop);
 
x_k(:,1) = b;
 
for k = 1:stop
    x_k(:,k+1) = x_k(:,k) + A_tilde_inv.*(b-A*x_k(:,k));
end





