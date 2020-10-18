%==============================================================================
% MSD Trapezoidal [ex4] (DOES NOT FULLY WORK YET)
% MOD13 4MC00 / Jelle Langedijk / TU/e
% TEMPLATE V1.1 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================
%% PREREQUISITIES
clear all; close all; clear vars; clc;

%% SYSTEM PARAMETERS
k =     10;
c =      2;
a1 =     0;
a2 =     1;

%% STEP SIZES
h1 = 2^(-1);
h2 = 2^(-2);
h3 = 2^(-3);

%% SYSTEM MATRIX
A = [0 1; -k -c];

%% MAKE q VECTOR
q(1,1) = [a1];          % q = [x,v]'
q(2,1) = [a2];          % q = [x,v]'
t = [0];

fac = 100;               % Accuracy Factor

for i = 1:10*(fac)
    mat = (1./(2*eye(2,2)-h3*A)).*(2*eye(2,2)+h3*A)*q(:,i);
    q(1,i+1) = mat(1,1);
    q(2,i+1) = mat(2,1);
    t(i+1) = i/fac;
end

%% PLOTS
figure(1)
plot(t,q(1,:));
 
 