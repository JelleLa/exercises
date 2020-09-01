%==============================================================================
% EX 3 MOD1 4MC00 / Jelle Langedijk / TU/e
% V1.0 (C) 2020 Jelle Langedijk, all rights reserved
% https://github.com/JelleLa/templates/tree/master/MATLAB
%==============================================================================

%% PREREQUISITIES
clear all; close all; clear vars;

%% DEFINITIONS

%NOTE= Du(x) = u'(x) hence du0 = u'(0)

%u(x) = e^(x)*cos(2x) +sin(sin(2x))
%Du(x) = e^x*cos(2x) - 2*e^x*sin(2x) + cos(sin(2x))*2cos(2x)
%Du(0) = 1 - 0 + 2 = 3
%Dhu(x) = u'(x) + O(h^2) --> O(h^2) = abs(Dhu(x) - Du(x))
%D~hu(x) = u'(x) + O(h) --> O(h) = abs(D~hu(x) - Du(x)) 

du0 = 3;
