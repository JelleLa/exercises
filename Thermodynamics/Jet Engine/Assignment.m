clear all;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation 
addpath("General"); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. Ready for use
global Runiv Pref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions. For the final assignment take the ones from the specific case you are supposed to do.                  
v1=200;Tamb=250;P3overP2=7;Pamb=55*kPa;mfurate=0.68*kg/s;AF=71.25;             % These are the ones from the book
cFuel='Gasoline';
%% Select all species
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);                                                          
Mi = [SpS.Mass];
%% Air comp
Xair = [0 0.21 0 0 0.79];                                                   % Order is important these are molefractions
MAir = Xair*Mi';
Yair = Xair.*Mi/MAir;
%% Fuel comp
Yfuel = Sp(iSp(1)).Elcomp;                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:3000];
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR 
    hia(:,i) = HNasa(TR,SpS(i)); %1e kolom gasooline etc.
    sia(:,i) = SNasa(TR,SpS(i));
end
hair_a= Yair*hia';                                                          % enthalpy of air for range of T 
sair_a= Yair*sia';                                                          % thermal part os entropy of air for range of T
whos hia sia hair_a sair_a                                                  % Shows dimensions of arrays on commandline
%% [1-2] Diffusor 
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;
T2 = interp1(hair_a,TR,h2);                                                 % Interpolated. How good is it?
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
s1thermal = Yair*si1';                                                      % Thermal part of entropy at T1
s2thermal = Yair*si2';                                                      % Thermal part of entropy at T2
lnPr = (s2thermal-s1thermal)/Rg;    
Pr = exp(lnPr);
P2 = P1*Pr;
s1  = s1thermal - Rg*log(P1/Pref);                                          % Specific Entropy at 1
s2  = s2thermal - Rg*log(P2/Pref);  

fprintf('\nStage %12s\n            | %9i %9i\n',sPart,1,2);
fprintf('%12s| %9.2f %9.2f\n','Temperature',T1,T2);
fprintf('%12s| %9.2f %9.2f  [kPa]\n','pressure',P1/kPa,P2/kPa);
fprintf('%12s| %9.2f %9.2f  [m/s]\n','velocity',v1,v2);
fprintf('%12s| %9.2f %9.2f  [J/kg]\n','enthalpy',h1,h2);
fprintf('%12s| %9.2f %9.2f  [J/kg]\n',' entropy',s1,s2)

%% [2-3] Compressor 

sPart = 'Compressor';

v3 = 0;
P3 = P3overP2 * P2;
s3 = s2;
s3thermal = s3 + Rg*log(P3/Pref);
T3 = interp1(sair_a,TR,s3thermal);

h3 = Yair * [0 ; HNasa(T3,Sp(iSp(2))) ; 0 ; 0 ; HNasa(T3,Sp(iSp(5))) ];

fprintf('\nStage %12s\n            | %9i %9i\n',sPart,2,3);
fprintf('%12s| %9.2f %9.2f\n','Temperature',T2,T3);
fprintf('%12s| %9.2f %9.2f  [kPa]\n','pressure',P2/kPa,P3/kPa);
fprintf('%12s| %9.2f %9.2f  [m/s]\n','velocity',v2,v3);
fprintf('%12s| %9.2f %9.2f  [J/kg]\n','enthalpy',h2,h3);
fprintf('%12s| %9.2f %9.2f  [J/kg]\n',' entropy',s2,s3)


%% [3-4] Combustor
sPart = 'Combustor';


v4=v3;
P4 = P3;

m_F_init = mfurate;
m_A_init = m_F_init * AF;
m_e = m_F_init + m_A_init;
m_tot_init = m_e;
m_O2_init = m_A_init * Yair(2);
m_N2_init = m_A_init * Yair(5);

Y_O2_init = (m_O2_init)/(m_tot_init);
Y_N2_init = (m_N2_init)/(m_tot_init);
Y_F_init = (m_F_init)/(m_tot_init);
Y_CO2_init = 0;
Y_H2O_init = 0;

n_F_init = (m_F_init)/(Sp(iSp(1)).Mass);
n_O2_init = (m_O2_init)/(Sp(iSp(2)).Mass); 
n_N2_init = (m_N2_init)/(Sp(iSp(5)).Mass); 
n_tot_init = n_N2_init + n_O2_init + n_F_init;

X_O2_init = (n_O2_init)/(n_tot_init);
X_N2_init = (n_N2_init)/(n_tot_init);
X_F_init = (n_F_init)/(n_tot_init);


n_F_final = 0;
n_O2_final = n_O2_init - n_F_init * 11.035;
n_N2_final = n_N2_init;
n_CO2_final = n_F_init*7.76;
n_H2O_final = n_F_init*6.55;
n_tot_final = n_F_final + n_O2_final + n_N2_final ...
    + n_CO2_final + n_H2O_final;

m_F_final = n_F_final*Sp(iSp(1)).Mass;
m_O2_final = n_O2_final*Sp(iSp(2)).Mass;
m_N2_final = n_N2_final*Sp(iSp(5)).Mass; 
m_CO2_final = n_CO2_final*Sp(iSp(3)).Mass;
m_H2O_final = n_H2O_final*Sp(iSp(4)).Mass;
m_tot_final = m_F_final + m_O2_final + m_N2_final ...
    + m_CO2_final + m_H2O_final;

Y_F_final = (m_F_final)/(m_tot_final);
Y_O2_final = (m_O2_final)/(m_tot_final);
Y_N2_final = (m_N2_final)/(m_tot_final);
Y_CO2_final = (m_CO2_final)/(m_tot_final);
Y_H2O_final = (m_H2O_final)/(m_tot_final);

X_F_final = (n_F_final)/(n_tot_final);
X_O2_final = (n_O2_final)/(n_tot_final);
X_N2_final = (n_N2_final)/(n_tot_final);
X_CO2_final = (n_CO2_final)/(n_tot_final);
X_H2O_final = (n_H2O_final)/(n_tot_final);


h_F = HNasa(T3,Sp(iSp(1)));
FA = AF^(-1);

h_RANGE = Y_F_final*HNasa(TR,Sp(iSp(1))) ...
+ Y_O2_final*HNasa(TR,Sp(iSp(2))) + Y_N2_final*HNasa(TR,Sp(iSp(5)))...
+ Y_CO2_final*HNasa(TR,Sp(iSp(3))) + Y_H2O_final*HNasa(TR,Sp(iSp(4)));

h4 = ((FA)/(1+FA))*h_F + ((1)/(1+FA))*h3; %see 8.2 page 20/564

T4 = interp1(h_RANGE, TR ,h4);

s4thermal = Y_F_final*SNasa(T4,Sp(iSp(1)))...
+ Y_O2_final*SNasa(T4,Sp(iSp(2))) + Y_N2_final*SNasa(T4,Sp(iSp(5)))...
+ Y_CO2_final*SNasa(T4,Sp(iSp(3))) + Y_H2O_final*SNasa(T4,Sp(iSp(4)));

Rg2 = (Runiv)/(X_O2_init*Sp(iSp(2)).Mass + X_N2_init*Sp(iSp(5)).Mass...
    + X_F_init*Sp(iSp(1)).Mass);
Rg3 = (Runiv)/(X_F_final*Sp(iSp(1)).Mass + X_O2_final*Sp(iSp(2)).Mass...
    + X_N2_final*Sp(iSp(5)).Mass + X_CO2_final*Sp(iSp(3)).Mass...
    + X_H2O_final*Sp(iSp(4)).Mass);

s4 = s4thermal - Rg3*log(P4/Pref);

ER = (((n_F_init)/(n_O2_init))/((1)/(11.035)))^(-1);

fprintf('\nStage %12s\n            | %9i %9i\n',sPart,3,4);
fprintf('%12s| %9.2f %9.2f\n','Temperature',T3,T4);
fprintf('%12s| %9.2f %9.2f  [kPa]\n','pressure',P3/kPa,P4/kPa);
fprintf('%12s| %9.2f %9.2f  [m/s]\n','velocity',v3,v4);
fprintf('%12s| %9.2f %9.2f  [J/kg]\n','enthalpy',h3,h4);
fprintf('%12s| %9.2f %9.2f  [J/kg]\n',' entropy',s3,s4);



%% [4-5] Turbine
sPart = 'Turbine';

s5=s4;

v5 = 0;

h5 = (m_A_init*(h2-h3) + m_tot_final*h4)/(m_tot_final);

T5 = interp1(h_RANGE, TR ,h5);

s5thermal = Y_F_final*SNasa(T5,Sp(iSp(1)))...
 + Y_O2_final*SNasa(T5,Sp(iSp(2))) + Y_N2_final*SNasa(T5,Sp(iSp(5)))...
 + Y_CO2_final*SNasa(T5,Sp(iSp(3))) + Y_H2O_final*SNasa(T5,Sp(iSp(4)));

P5  = Pref*exp((s5-s5thermal)/(-Rg3));

fprintf('\nStage %12s\n            | %9i %9i\n',sPart,4,5);
fprintf('%12s| %9.2f %9.2f\n','Temperature',T4,T5);
fprintf('%12s| %9.2f %9.2f  [kPa]\n','pressure',P4/kPa,P5/kPa);
fprintf('%12s| %9.2f %9.2f  [m/s]\n','velocity',v4,v5);
fprintf('%12s| %9.2f %9.2f  [J/kg]\n','enthalpy',h4,h5);
fprintf('%12s| %9.2f %9.2f  [J/kg]\n',' entropy',s4,s5);

%% [5-6] Nozzle
sPart = 'Nozzle';

v5=v4;
P6 =Pamb;
s6 = s5;
s6thermal = s6 + Rg3*log(P6/Pref);

s6_range = Y_F_final*SNasa(TR,Sp(iSp(1)))...
+ Y_O2_final*SNasa(TR,Sp(iSp(2))) + Y_N2_final*SNasa(TR,Sp(iSp(5)))...
+ Y_CO2_final*SNasa(TR,Sp(iSp(3))) + Y_H2O_final*SNasa(TR,Sp(iSp(4)));

T6 = interp1(s6_range,TR,s6thermal);

h6 = Y_F_final*HNasa(T6,Sp(iSp(1))) + Y_O2_final*HNasa(T6,Sp(iSp(2)))...
+ Y_N2_final*HNasa(T6,Sp(iSp(5))) + Y_CO2_final*HNasa(T6,Sp(iSp(3)))...
+ Y_H2O_final*HNasa(T6,Sp(iSp(4)));

v6 = sqrt(2*(h5-h6));


fprintf('\nStage %12s\n            | %9i %9i\n',sPart,5,6);
fprintf('%12s| %9.2f %9.2f\n','Temperature',T5,T6);
fprintf('%12s| %9.2f %9.2f  [kPa]\n','pressure',P5/kPa,P6/kPa);
fprintf('%12s| %9.2f %9.2f  [m/s]\n','velocity',v5,v6);
fprintf('%12s| %9.2f %9.2f  [J/kg]\n','enthalpy',h5,h6);
fprintf('%12s| %9.2f %9.2f  [J/kg]\n',' entropy',s5,s6);



