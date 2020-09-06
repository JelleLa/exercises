%==============================================================================
% MOD2 4MC00 / Jelle Langedijk / TU/e
% V1.0 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================

%% PREREQUISITIES
clear all; close all; clear vars;

%% MATRIX FILLER

N = 8;


Ah = zeros(N);

for i = 1:N
    Ah(i,i) = -2;
end   

for j = 1:N-1
    Ah(j,j+1) = +1;
end  

for k = 2:N
    Ah(k,k-1) = +1;
end  

%% 

