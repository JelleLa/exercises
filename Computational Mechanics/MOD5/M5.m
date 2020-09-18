%==============================================================================
% MOD4 4MC00 / Jelle Langedijk / TU/e
% TEMPLATE V1.1 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================

%% PREREQUISITIES
clear all; close all; clear vars; clc;

%% EX 3
L = 4;

x = [0:(L/50):L];
x1 = 0*L;
x2 = 1/3*L;
x3 = 2/3*L;
x4 = L;

u = [1 3 3 2];

N1 = (1-(3*x)/(L)).*(1-(3*x)/(2*L)).*(1-(x)/(L));
N2 = (1-(3*x)/(2*L)).*(1-(x)/(L)).*x.*((9)/(L));
N3 = (1-(x)/(L)).*(1-(3*x)/(L)).*x.*((-9)/(2*L));
N4 = (1-(3*x)/(L)).*(1-(3*x)/(2*L)).*x.*((1)/(L));

func = N1*u(1) + N2*u(2) + N3*u(3) + N4*u(4);

figure(1)
plot(x, N1);
hold on
plot(x, N2);
hold on
plot(x, N3);
hold on
plot(x, N4);
hold on
plot(x, func);
hold on
plot(x1, u(1), "r*");
hold on
plot(x2, u(2), "r*");
hold on
plot(x3, u(3), "r*");
hold on
plot(x4, u(4), "r*");
hold on
legend("N1","N2","N3","N4","func")

%% EX 4

% n = p+1 , p = n-1

clear all; close all; clear vars; clc;

DOF = 20;
tab = zeros(DOF,2);
tab2 = zeros(DOF,2);

for n = 1:DOF

[E,c] = polyapprox(n-1);
tab(n,1) = n;
tab(n,2) = E;

end

for n = 1:DOF

[E,c] = polyapprox(n-1);
tab2(n,1) = n;
tab2(n,2) = c;

end

figure(2)
semilogy(tab(:,1), tab(:,2))
title("E versus n")

figure(3)
semilogy(tab2(:,1), tab2(:,2))
title("c versus n")



