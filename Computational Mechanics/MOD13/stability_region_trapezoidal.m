%==============================================================================
% STABILITY REGION TRAPEZOIDAL [ex3]
% MOD13 4MC00 / Jelle Langedijk / TU/e
% TEMPLATE V1.1 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================
%% PREREQUISITIES
clear all; close all; clear vars; clc;

%% Define Imaginary Unit
i=sqrt(-1);

%% System Parameters
k = 10;
c = 2;

%% Calculate Eigenvalues

labda1 = (-c + sqrt(c^2-4*k))/(2);
labda2 = (-c - sqrt(c^2-4*k))/(2);

%% Arrays of labda*h

ev1 = [labda1*2^(-1) labda1*2^(-2) labda1*2^(-3)];
ev2 = [labda2*2^(-1) labda2*2^(-2) labda2*2^(-3)];


%% Specify real range and number of points
lhr0=-3;
lhr1=3;
Nr = 1024;

%% Specify imaginary range and number of points
lhi0=-3;
lhi1=3;
Ni = 1024;

%% Construct mesh
lhr = linspace(lhr0,lhr1,Nr);
lhi = linspace(lhi0,lhi1,Ni);
[x,y] = meshgrid(lhr,lhr);

%% Calculate lh
lh=x+i*y;

%% Amplification factor of Trapezoid
theta = (1./((1-(lh/2))).*(1+(lh/2)));

%% Plot contours of magnitude(theta)=1
contour(x,y,abs(theta),[0:0.05:1],'b-');
hold on;

plot(real(ev1(1)), imag(ev1(1)), "r*");
text(real(ev1(1)), imag(ev1(1)), "h = 2^{-1}");
hold on;
plot(real(ev1(2)), imag(ev1(2)), "r*");
text(real(ev1(2)), imag(ev1(2)), "h = 2^{-2}");
hold on;
plot(real(ev1(3)), imag(ev1(3)), "r*");
text(real(ev1(3)), imag(ev1(3)), "h = 2^{-3}");
hold on;
plot(real(ev2(1)), imag(ev2(1)), "r*");
text(real(ev2(1)), imag(ev2(1)), "h = 2^{-1}");
hold on;
plot(real(ev2(2)), imag(ev2(2)), "r*");
text(real(ev2(2)), imag(ev2(2)), "h = 2^{-2}");
hold on;
plot(real(ev2(3)), imag(ev2(3)), "r*");
text(real(ev2(3)), imag(ev2(3)), "h = 2^{-3}");


axis([lhr0,lhr1,lhi0,lhi1]);
axis('square');
xlabel('\lambda h_r');
ylabel('\lambda h_i');
grid on