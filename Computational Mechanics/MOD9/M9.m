%==============================================================================
% MOD9 4MC00 / Jelle Langedijk / TU/e
% TEMPLATE V1.1 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================



%% EX 4
% PREREQUISITIES
clear all; close all; clear vars; clc;

% DEFINE CARTESIAN BASE
e = cartesianbasis2d('e1', 'e2');
e1 = e(1);
e2 = e(2);

% PRODUCT TEST
prod1 = dot(e1, e1);
prod2 = dot(e1, e2);

% DEFINE x AND u
x = 4*e1 + 2*e2;
u = -1*e1 + 2*e2;

prod3 = dot(x, e1);
prod4 = dot(x, e2);

prod5 = dot(x, u);

N  = norm(u);

% QUIVER PLOT
quiver(x, u, 0)

% DIADIC PRODUCT
diad1 = x * u; 

%% EX 5
% PREREQUISITIES
clear all; close all; clear vars; clc;

e = cartesianbasis2d('e1', 'e2');
e1 = e(1);
e2 = e(2);

a = [3*e1 - 1*e2 ; e1 + 3*e2 ; 5*e1 ];

C = 2*e1*e1 - e1*e2 -e2*e1 + e2*e2;

at = a' ;

prod1 = dot(C, a) ;
prod2 = dot(at, prod1) ;

O = 0*e1*e1 + 0*e1*e2 + 0*e2*e1 + 0*e2*e2;

D = [C O O; O -1*C O ; O O (1/5)*C];

prod3 = dot(D, a);



