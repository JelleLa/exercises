%==============================================================================
% MOD2 4MC00 / Jelle Langedijk / TU/e
% V1.1 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================

%% PREREQUISITIES
clear all; close all; clear vars; clc;


 

%% VERTICAL LOAD
% f  = 10^(2)*e^(-1000*(x-(1)/(2))^(2))

h = 2^(-3);             %Step Size
L = 1;                  %Length of xi (See BC's: from 0 to 1)
xi = [0:h:L];
fi = 10.^(2)*exp(-1000*(xi-(1)/(2)).^(2));
b = fi';

%% MATRIX FILLER

N = length(xi) ;


A = zeros(N);

for i = 1:N
    A(i,i) = -2;
end   

for j = 1:N-1
    A(j,j+1) = +1;
end  

for k = 2:N
    A(k,k-1) = +1;
end 
    
%% SOLVE
Ah = (-1)/(h^(2)).*A ;  %Ah Matrix

U = Ah\b;

%% MIDPOINT DISPLACEMENT

mid = ceil(length(xi)/2);
disp = b(mid) - U(mid);



%% PLOTS

figure(1)
plot(xi , U) ;
xlabel("x_i") ;
ylabel("u^h_i") ;
title("MODULE 2 EX 3a") ;

figure(2)
plot(xi , fi) ;
xlabel("x_i") ;
ylabel("f_i") ;
title("MODULE 2 EX 3b") ;


