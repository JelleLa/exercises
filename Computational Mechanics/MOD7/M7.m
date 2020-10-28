%==============================================================================
% MOD7 4MC00 / Jelle Langedijk / TU/e
% TEMPLATE V1.1 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================

%% EX 3A
%Kgg = 1
%
clear all; close all; clear vars; clc;

K = [   1, -1, 0, 0, 0; 
        -1, 2, -1, 0, 0;
        0, -1, 2, -1, 0;
        0, 0, -1 , 2, -1;
        0, 0, 0, -1, 1;
        ];

n = length(K);    
Kgg = K(1  , 1  );
Kfg = K(2:n, 1  );
Kgf = K(1  , 2:n);
Kff = K(2:n, 2:n);

ug  = 0;
qf  = transpose([0 0 0 1]);

% solve the system for unknowns
uf = Kff \ (qf - Kfg*ug);
u  = [ug; uf];

qg = Kgg * ug + Kgf * uf;
q = [qg; qf];    
    
    
%% EX 3B
clear all; close all; clear vars; clc;
K = [   1, -1, 0, 0, 0; 
        -1, 2, -1, 0, 0;
        0, -1, 2, -1, 0;
        0, 0, -1 , 2, -1;
        0, 0, 0, -1, 1;
        ];
n = length(K);    
Kgg = K(1  , 1  );
Kgf = K(1, 2:4);
Kgh = K(1, 5);
Kfg = (K(2:4, 1));
Kff = K(2:4, 2:4);
Kfh = (K(2:4, 5));
Khg = K(5,1);
Khf = K(5, 2:4);
Khh = K(5,5);

ug = 0;
uh = 1;
qf = [0;0;0];


uf = (Kff)\(qf - Kfg*ug - Kfh*uh) ;
u = [ug; uf; uh]; 

qg = Kgg*ug + Kgf*uf + Kgh*uh;
qh = Khg*ug + Khf*uf + Khh*uh;

q = [qg;qf;qh];

%% EX 3C
clear all; close all; clear vars; clc;
K = [   1, -1, 0, 0, 0; 
        -1, 2, -1, 0, 0;
        0, -1, 2, -1, 0;
        0, 0, -1 , 2, -1;
        0, 0, 0, -1, 1;
        ];
n = length(K); 

Kgg = K(1:(n-1),1:(n-1));
Kgf = K(1:n-1 , n);
Kfg = K(n ,1:n-1);
Kff = K(n,n);

qg = transpose([0 0 0 0]) ;
uf = 1;


ug = Kgg\(qg - Kgf*uf) ;

qf = Kfg*ug +Kff*uf;

u = [ug' , uf]';
q = [qg' , qf]';


    
    


