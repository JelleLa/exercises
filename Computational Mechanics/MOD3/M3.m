%==============================================================================
% MOD3 4MC00 / Jelle Langedijk / TU/e
% V1.1 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================

%% PREREQUISITIES
clear all; close all; clear vars; clc;


 

%% VERTICAL LOAD
N = 9;
h = 1/(N);             
L = 1;                  
xi = [0:h:L];


%% MATRIX FILLER (Neumann)
A = zeros(N);
A(1,1) = -(3)/(2*h);
A(1,2) = (2)/(h);
A(1,3) = -(1)/(2*h);

A(N,N-2) = (1)/(2*h);
A(N,N-1) = -(2)/(h);
A(N,N) = (3)/(2*h);

for i = 2:N-1
    A(i,i-1) = -(1)/(h^(2));
end   

for j = 2:N-1
    A(j,j) = (2)/(h^(2));
end  

for k = 2:N-1
    A(k,k+1) = -(1)/(h^(2));
end 
    
%% SOLVE

for i = 1:N
    b(i) = 2;
end 

b = transpose(b);
U = A\b;

v = null(A);
%% MATRIX FILLER (Dirichlet)
for g = 1:N
    A(g,1) = 0;
end

UD = A\b;
