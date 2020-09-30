%==============================================================================
% MOD7 4MC00 / Jelle Langedijk / TU/e
% TEMPLATE V1.1 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================

%% PREREQUISITIES
clear all; close all; clear vars; clc;

%% ASSEMBLY

labda = 80;
R0 = 20;
Ri = 2;
m = 5;
h = (R0-Ri)/m;
Ti = 360;
T0 = 300;

K = zeros((R0-Ri)/h,(R0-Ri)/h);

re0 = Ri;
re1 = Ri +h;

for i = 1:m %(R0-Ri)/h

Ke = 0.5*labda*(re1+re0)/(re1-re0).*[1 -1 ; -1 1];

K(i, i) = Ke(1,1);
K(i, i+1) = Ke(1,2);
K(i+1, i) = Ke(2,1);
K(i+1, i+1) = sum([Ke(2,2) K(i,i)]);  %WAAROM TELT DIT NIET OP VOOR IEDERE ITERATIE?

re0 = re0+h;
re1 = re1+h;
end

%% PARTITIONING 
% Qi and Q0 are my unknowns
% Ti and T0 given --> partition problem identical to M7 3b. 
% u = T

n = length(K);   

Kgg = K(1,1);
Kgf = K(1,2:(n-1));
Kgh = K(1,n);
Kfg = (K(2:(n-1), 1));
Kff = K(2:(n-1), 2:(n-1));
Kfh = (K(2:(n-1), n));
Khg = K(n,1);
Khf = K(n, 2:(n-1));
Khh = K(n,n);

ug = Ti;
uh = T0;
qf = zeros((n-2),1);


uf = (Kff)\(qf - Kfg*ug - Kfh*uh) ;
u = [ug; uf; uh]; 

qg = Kgg*ug + Kgf*uf + Kgh*uh;
qh = Khg*ug + Khf*uf + Khh*uh;

q = [qg;qf;qh];





    
    


