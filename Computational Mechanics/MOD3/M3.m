%==============================================================================
% MOD3 4MC00 / Jelle Langedijk / TU/e
% V1.1 (C) 2020 Jelle Langedijk, all rights reserved
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

%% RMS (expression Lu determined by hand). Also see slides lecture 2. Does not work yet

% u_xi = (1/2)-(1/2)*cos(2*pi*xi);
% Lu_xi = -2*(pi)^(2)*cos(2*pi*xi);
% 
% RMS = rms(transpose(rms(transpose(Lu_xi)- A.*transpose(u_xi))));



%% SOLVE

for i = 1:N
    b(i) = 2;
end 

b = transpose(b);
U = A\b;

v = null(A);
%% MATRIX FILLER (Dirichlet)
A(1,1) = 0;

UD = A\b;
