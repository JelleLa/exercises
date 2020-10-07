%==============================================================================
% MOD10 4MC00 / Jelle Langedijk / TU/e
% TEMPLATE V1.1 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================

%% EX 1

% PREREQUISITIES
clear all; close all; clear vars; clc;

e  = cartesianbasis2d('e1', 'e2');
e1 = e(1);
e2 = e(2);
xie = [ -e1-e2 e1-e2 e1+e2 -e1+e2 ];

Q   = [ 1  2  3  4 ];

 %femplot(xie, Q, 'Nodes', 'on')

xi_0  = 0*e1 + 0*e2;        % Origin
xi_P  = (1/2)*e1 + 1*e2 ;   % P
xi_1 = -1*e1 + -1*e2;       % Node 1
xi_2 = 1* e1 -1*e2;         % Node 2
xi_3 = 1* e1 + 1* e2;       % Node 3
xi_4 = -1* e1 + 1* e2;      % Node 4

xi = xi_4;
xi1 = dot(xi, e1);
xi2 = dot(xi, e2);
Ne  = [ 1/4*(1-xi1)*(1-xi2)
        1/4*(1+xi1)*(1-xi2)
        1/4*(1+xi1)*(1+xi2)
        1/4*(1-xi1)*(1+xi2) ];
    
grad_Ne  = [ -1/4*(1-xi2)*e1 -1/4*(1-xi1)*e2
             1/4*(1-xi2)*e1 -1/4*(1+xi1)*e2
             1/4*(1+xi2)*e1 +1/4*(1+xi1)*e2
             -1/4*(1+xi2)*e1 +1/4*(1-xi1)*e2 ];
         
%% EX 2 (Lagrange)

% PREREQUISITIES
clear all; close all; clear vars; clc;

e  = cartesianbasis2d('e1', 'e2');
e1 = e(1);
e2 = e(2);
xie = [ -e1-e2 e1-e2  e1+e2 -e1+e2 -e2  e1 +e2 -1*e1 0*e1+0*e2 ];

QL   = [ 1  2  3  4  5  6  7  8  9 ];

femplot(xie, QL, 'Nodes', 'on')

xi_0  = 0*e1 + 0*e2;        % Origin
xi_P  = (1/2)*e1 + 1*e2 ;   % P

% See Handwritten derivation and Dictation
xi = xi_P;
xi1 = dot(xi, e1);
xi2 = dot(xi, e2);
Ne  = [ 1/4*(1-xi1)*abs(xi1)*(1-xi2)*abs(xi2)
        -1/4*(1+xi1)*abs(xi1)*(1-xi2)*abs(xi2)
        1/4*(1+xi1)*abs(xi1)*(1+xi2)*abs(xi2)
        -1/4*(1-xi1)*abs(xi1)*(1+xi2)*abs(xi2)
        -1/2*(1-xi1)*(1+xi1)*(1-xi2)*abs(xi2)
        1/2*(1-xi2)*(1+xi1)*(1+xi2)*abs(xi1)
        1/2*(1+xi2)*(1+xi1)*(1-xi1)*abs(xi2)
        -1/2*(1+xi2)*(1-xi2)*(1-xi1)*abs(xi1)
        (1+xi2)*(1-xi2)*(1+xi1)*(1-xi1) ];

u_e = [11 14 12 17 15 13 18 16 19]';

u_h = transpose(Ne)*u_e;

%% EX 5a (Weird Element)

% PREREQUISITIES
clear all; close all; clear vars; clc;

e  = cartesianbasis2d('e1', 'e2');
e1 = e(1);
e2 = e(2);
xie = [ (-e1-e2); (e1-e2);  (e1+e2); (-e1+e2); (-(1/3)*e1-e2); (+(1/3)*e1-e2); (e1-(1/3)*e2); (e1+(1/3)*e2); (+(1/3)*e1+e2); (-(1/3)*e1+e2); (-e1+(1/3)*e2); (-e1-(1/3)*e2); (-(1/3)*e1-(1/3)*e2); (+(1/3)*e1-(1/3)*e2); (+(1/3)*e1+(1/3)*e2); (-(1/3)*e1+(1/3)*e2);];

QL   = [ 1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16];

%femplot(xie, QL, 'Nodes', 'on')

%% EX 5e (Weird Element)
% PREREQUISITIES
clear all; close all; clear vars; clc;

xi1 = [-1:0.001:1];
xi2 = 1;

% See handwritten derivations of N4 and N16
N4 = -(81/256).*(xi1+1/3).*(xi2 + 1/3).*(xi1-1/3).*(xi2-1/3).*(xi1-1).*(xi2+1);

N16 = (729/256).*(xi1-1).*(xi2-1).*(xi1-1/3).*(xi2+ 1/3).*(xi1+1).*(xi2+1);

figure(1)
plot(xi1, N4)
hold on
plot(xi1, N16)
xlabel("$\xi_{1}$", 'Interpreter', 'latex')
ylabel("$N$", 'Interpreter', 'latex')
title("$\xi_{2} = 1$", 'Interpreter', 'latex')
legend("N1", "N16")

%% EX 5f (Weird Element)
% PREREQUISITIES
clear all; clear vars; clc; 

xi1 = [-1:0.001:1];
xi2 = 1/3;

% See handwritten derivations of N4 and N16
N4 = -(81/256).*(xi1+1/3).*(xi2 + 1/3).*(xi1-1/3).*(xi2-1/3).*(xi1-1).*(xi2+1);

N16 = (729/256).*(xi1-1).*(xi2-1).*(xi1-1/3).*(xi2+ 1/3).*(xi1+1).*(xi2+1);

figure(2)
plot(xi1, N4)
hold on
plot(xi1, N16)
xlabel("$\xi_{1}$", 'Interpreter', 'latex')
ylabel("$N$", 'Interpreter', 'latex')
title("$\xi_{2} = \frac{1}{3}$", 'Interpreter', 'latex')
legend("N1", "N16")
