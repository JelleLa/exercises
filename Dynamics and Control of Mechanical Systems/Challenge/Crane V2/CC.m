% -------------------------------------------------------------------------
% 4DB00 Dynamics and control of mechanical systems 2019-2020
% Crane challenge.
% -------------------------------------------------------------------------
clc;close all;clearvars
addpath Toolbox 
addpath Toolbox/Functions

%% Symbols
% Define the degrees of freedom and their first two derivative w.r.t. time
% as symbolic expressions
syms phi1 phi2 x1 dphi1 dphi2 dx1 ddphi1 ddphi2 ddx1   

% Define the system parameters, actuator force and time 
% as symbolic expressions
syms g L0 L1 L2 L m0 m1 m2 m3 kh dW dh    
syms FA                                                
syms t                  

% Place degrees of freedom (and their derivatives) in the generalized
% coordinates and the prescribed displacements.
Components.q        = [x1; phi2];
Components.dq       = [dx1; dphi2];
Components.ddq      = [ddx1; ddphi2];
Components.s_pd     = [phi1];
Components.ds_pd    = [dphi1];
Components.dds_pd   = [ddphi1];
Components.F        = [FA];
Components.time     = [t];

% Define system parameters as symbols
Components.parameters.g.symbol  = g;
Components.parameters.L0.symbol = L0;
Components.parameters.L1.symbol = L1;
Components.parameters.L2.symbol = L2;
Components.parameters.m0.symbol = m0;
Components.parameters.m1.symbol = m1;
Components.parameters.m2.symbol = m2;
Components.parameters.m3.symbol = m3;
Components.parameters.kh.symbol = kh;
Components.parameters.dW.symbol = dW;
Components.parameters.dh.symbol = dh;

%% ------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------> Please complete the following code for Deliverable 1:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

Components.parameters.g.value  = [9.81];    %m/(s^2)
Components.parameters.L0.value = [20];    %m
Components.parameters.L1.value = [25];    %m
Components.parameters.L2.value = [10];    %m
Components.parameters.m0.value = [10000];    %kg
Components.parameters.m1.value = [10000];    %kg
Components.parameters.m2.value = [2000];    %kg
Components.parameters.m3.value = [500];    %kg
Components.parameters.kh.value = [0];    %N/m
Components.parameters.dW.value = [100];    %N*s/(rad*m) 200
Components.parameters.dh.value = [200];    %N*s/m

% Insert symbolic expressions for the kinetic and potential energy 
% and the generalized forces:
T   = [0.5*dx1^2*(m2 + m3) + 0.5*x1^2*dphi1^2*(m2 + m3) + 0.5*L2^2*dphi2^2*m3 + m3*dx1*L2*dphi2*(cos(phi1)*sin(phi2) - cos(phi2)*sin(phi1)) - m3*x1*L2*dphi1*dphi2*(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)) + (1/6)*m1*L1^2*dphi1^2];   % Kinetic energy
V   = [m0*g*0.5*L0 + m1*g*(L0+0.5*L1*sin(phi1)) + m2*g*(L0+ x1*sin(phi1)) + m3*g*(L0 + x1*sin(phi1) - L2*sin(phi2)) + 0.5*kh*(x1-0.5*L1)^2];   % Potential energy
Qnc = [FA - dh*dx1 ; -dW*L2^2*dphi2];   % Generalised forces

% Calculate the non-linear equations of motion (statement is complete):
Nonlinear_EOM = DOMS('NonlinearEOM',Components,T,V,Qnc);

% Insert the settings for the simulation:
t_sim     = [100];  % Simulation time
q_init    = [11.5 ; pi/3];  % Initial value for q
dq_init   = [0;0];  % Initial value for dq
s_pd_of_t = [0];  % Prescribed displacement a function of time t (e.g. sin(t))
FA_of_t   = [0];  % External force as a function of time t (e.g. cos(t))

% Simulate the non-linear equations of motion (statement is complete):
Sim_Nonlinear = DOMS('Sim_Nonlinear',Components, Nonlinear_EOM,...
                                   t_sim,q_init,dq_init,s_pd_of_t,FA_of_t);

% Simulate the virtual test setup in Simulink (statement is complete):
Sim_Setup  = DOMS('Sim_Setup',Components,t_sim,...
                                   q_init,dq_init,s_pd_of_t,FA_of_t);

%--------------------------------------------------------------------------
%--------> Please generate the plots for questions e and f here:
%--------------------------------------------------------------------------

str1 = "$$ Time [s] $$";
str2 = "$$ x1 [m] $$";
str3 = "$$ \varphi_{2} [rad] $$";


figure()
hold on
plot(Sim_Nonlinear.t,Sim_Nonlinear.qt(:,1))
plot(Sim_Setup.t,Sim_Setup.qt(:,1))
xlim([0 100])
ylim([9 14])
title("Position of the hoist block" , 'Interpreter','latex')
xlabel(str1,'Interpreter','latex')
ylabel(str2,'Interpreter','latex')
legend('Derived Model' , 'Simulation')


figure()
hold on
plot(Sim_Nonlinear.t,Sim_Nonlinear.qt(:,2))
plot(Sim_Setup.t,Sim_Setup.qt(:,2))
xlim([0 100])
ylim([1 2])
title("Absolute angle of the chain", 'Interpreter','latex')
xlabel(str1,'Interpreter','latex')
ylabel(str3,'Interpreter','latex')
legend('Derived Model' , 'Simulation')

return

%% ------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------> Please complete the following code for Deliverable 2:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Insert the linearization points
q_0_stable   = [0.5*25 ; (pi/2)];    %Stable equilibrium for linearization


% Linearize the equations of motion (statements are complete):
LinearEOM_stable   = DOMS('Linearize',Components,Nonlinear_EOM,q_0_stable);  


% Insert the settings for the simulation:
t_sim             = [100];   % Simulation time
q_init            = [11.5 ; pi/3];
dq_init           = [0 ; 0];
q1_init_stable    = [q_init - q_0_stable];   % Initial value for q1
dq1_init_stable   = [0;0];   % Initial value for dq1
FA_of_t           = [0];   % External force as a function of time t (e.g. cos(t))


% Simulate the linear equations of motion (statements are complete):
Sim_Linear_stable   = DOMS('Sim_Linear',Components,LinearEOM_stable,...
                            t_sim,q1_init_stable,dq1_init_stable, FA_of_t);



%--------------------------------------------------------------------------
%--------> Please generate the plots for questions d and e here:
%--------------------------------------------------------------------------

%%NOT COMPARED
%Stable Plots
    %x1
    figure()
    hold on
    plot(Sim_Linear_stable.t, Sim_Linear_stable.q1t(:,1) + q_0_stable(1) )
    xlim([0 100])
    ylim([9 14])
    title("Position of the hoist block (Stable)" , 'Interpreter','latex')
    xlabel(str1,'Interpreter','latex')
    ylabel(str2,'Interpreter','latex')
    legend('Derived Model' , 'Simulation')
    
    %phi2
    figure()
    hold on
    plot(Sim_Linear_stable.t, Sim_Linear_stable.q1t(:,2) + q_0_stable(2) )
    xlim([0 100])
    ylim([1 2])
    title("Absolute angle of the chain (Stable)", 'Interpreter','latex')
    xlabel(str1,'Interpreter','latex')
    ylabel(str3,'Interpreter','latex')
    legend('Derived Model' , 'Simulation')
    
%Unstable Plots
    %x1
    figure()
    hold on
    plot(Sim_Linear_unstable.t, Sim_Linear_unstable.q1t(:,1) + q_0_unstable(1) )
    xlim([0 100])
    ylim([9 14])
    title("Position of the hoist block (Unstable)" , 'Interpreter','latex')
    xlabel(str1,'Interpreter','latex')
    ylabel(str2,'Interpreter','latex')
    legend('Derived Model' , 'Simulation')
    
    %phi2
    figure()
    hold on
    plot(Sim_Linear_unstable.t, Sim_Linear_unstable.q1t(:,2) + q_0_unstable(2) )
    xlim([0 100])
    ylim([1 2])
    title("Absolute angle of the chain (Unstable)", 'Interpreter','latex')
    xlabel(str1,'Interpreter','latex')
    ylabel(str3,'Interpreter','latex')
    legend('Derived Model' , 'Simulation')

%%COMPARED
%Stable Plots
    %x1
    figure()
    hold on
    plot(Sim_Linear_stable.t, Sim_Linear_stable.q1t(:,1) + q_0_stable(1) )
    plot(Sim_Setup.t,Sim_Setup.qt(:,1))
    xlim([0 100])
    ylim([9 14])
    title("Position of the hoist block (Stable)" , 'Interpreter','latex')
    xlabel(str1,'Interpreter','latex')
    ylabel(str2,'Interpreter','latex')
    legend('Derived Model' , 'Simulation')
    
    %phi2
    figure()
    hold on
    plot(Sim_Linear_stable.t, Sim_Linear_stable.q1t(:,2) + q_0_stable(2) )
    plot(Sim_Setup.t,Sim_Setup.qt(:,2))
    xlim([0 100])
    ylim([1 2])
    title("Absolute angle of the chain (Stable)", 'Interpreter','latex')
    xlabel(str1,'Interpreter','latex')
    ylabel(str3,'Interpreter','latex')
    legend('Derived Model' , 'Simulation')
    
%Unstable Plots
    %x1
    figure()
    hold on
    plot(Sim_Linear_unstable.t, Sim_Linear_unstable.q1t(:,1) + q_0_unstable(1) )
    plot(Sim_Setup.t,Sim_Setup.qt(:,1))
    xlim([0 100])
    ylim([9 14])
    title("Position of the hoist block (Unstable)" , 'Interpreter','latex')
    xlabel(str1,'Interpreter','latex')
    ylabel(str2,'Interpreter','latex')
    legend('Derived Model' , 'Simulation')
    
    %phi2
    figure()
    hold on
    plot(Sim_Linear_unstable.t, Sim_Linear_unstable.q1t(:,2) + q_0_unstable(2) )
    plot(Sim_Setup.t,Sim_Setup.qt(:,2))
    xlim([0 100])
    ylim([1 2])
    title("Absolute angle of the chain (Unstable)", 'Interpreter','latex')
    xlabel(str1,'Interpreter','latex')
    ylabel(str3,'Interpreter','latex')
    legend('Derived Model' , 'Simulation')
    
%%D3-------------------------------------------------------------------
    
w=sqrt(eig(LinearEOM_stable.K , LinearEOM_stable.M));
w1 = w(1);
w2= w(2);

% Insert the settings for the simulation:
t_sim     = [200];  % Simulation time
q_init    = [12.5 ; pi/2];  % Initial value for q
dq_init   = [0;0];  % Initial value for dq
s_pd_of_t = [0];  % Prescribed displacement a function of time t (e.g. sin(t))
FA_of_t   = [1500*sin(1.3*w2*t)];  % VARIATE PER QUESTION

% Simulate the virtual test setup in Simulink (statement is complete):
Sim_Setup  = DOMS('Sim_Setup',Components,t_sim,...
                                   q_init,dq_init,s_pd_of_t,FA_of_t);

figure()
hold on
plot(Sim_Setup.t,Sim_Setup.qt(:,1))
xlim([0 200])
ylim([10 25])
title("Position of the hoist block" , 'Interpreter','latex')
xlabel(str1,'Interpreter','latex')
ylabel(str2,'Interpreter','latex')

figure()
hold on
plot(Sim_Setup.t,Sim_Setup.qt(:,2))
xlim([0 200])
ylim([1 2])
title("Absolute angle of the chain", 'Interpreter','latex')
xlabel(str1,'Interpreter','latex')
ylabel(str3,'Interpreter','latex')


G = DOMS("Make_tf" , LinearEOM_stable);

figure()
bode(G);

E = abs(evalfr(G, i*w2));

Gl = G(1) + Components.parameters.L2.value * G(2)

%%D4-------------------------------------------------------------------
figure()
bode(Gl);

Z = zero(Gl);
P = pole(Gl);

figure()
step(40*Gl);
xlim([0 100]);




% Insert the settings for the simulation:
t_sim     = [100];  % Simulation time
q_init    = [12.5 ; pi/2];  % Initial value for q
dq_init   = [0;0];  % Initial value for dq
s_pd_of_t = [0];  % Prescribed displacement a function of time t (e.g. sin(t))
FA_of_t   = [40];  % VARIATE PER QUESTION

% Simulate the virtual test setup in Simulink (statement is complete):
Sim_Setup  = DOMS('Sim_Setup',Components,t_sim,...
                                   q_init,dq_init,s_pd_of_t,FA_of_t);

str4 = "$$ x_{l} [m] $$";       
                               
figure()
plot(Sim_Setup.t, Sim_Setup.qt(:,1) - Components.parameters.L2.value*cos(Sim_Setup.qt(:,2)))
xlim([0 100]);
ylim([10 25]);
title("Position of the load" , 'Interpreter','latex')
xlabel(str1,'Interpreter','latex')
ylabel(str4,'Interpreter','latex')

C = 40;

HA = (C*Gl)/(1+ C*Gl);

HA_min = minreal(HA);

figure()
step(HA);
xlim([0 300]);



t_sim     = [300];  % Simulation time
q_init    = [12.5 ; pi/2];  % Initial value for q
dq_init   = [0;0];  % Initial value for dq
s_pd_of_t = [0];  % Prescribed displacement a function of time t (e.g. sin(t))
P1 = 4;
P2 = 40;
P3 = 400;
reference_signal = 'Step';
reference_settings = [12.5 , 13.5 , 0];
Sim_step_P1 = DOMS("Sim_Control", Components, t_sim, q_init , dq_init , s_pd_of_t, P1 , reference_signal , reference_settings);
Sim_step_P2 = DOMS("Sim_Control", Components, t_sim, q_init , dq_init , s_pd_of_t, P2 , reference_signal , reference_settings);
Sim_step_P3 = DOMS("Sim_Control", Components, t_sim, q_init , dq_init , s_pd_of_t, P3 , reference_signal , reference_settings);


figure()
plot(Sim_step_P1.t, Sim_step_P1.qt(:,1) - Components.parameters.L2.value*cos(Sim_step_P1.qt(:,2)));
hold on
plot(Sim_step_P2.t, Sim_step_P2.qt(:,1) - Components.parameters.L2.value*cos(Sim_step_P2.qt(:,2)));
hold on
plot(Sim_step_P3.t, Sim_step_P3.qt(:,1) - Components.parameters.L2.value*cos(Sim_step_P3.qt(:,2)));
xlim([0 300]);
ylim([11 15]);
title("Position of the load" , 'Interpreter','latex');
xlabel(str1,'Interpreter','latex');
ylabel(str4,'Interpreter','latex');
legend('P=4' , 'P=40' , 'P=400');

%%D5-------------------------------------------------------------------

%Lead-Compensators

s = tf("s");
C1 = 61.944*((8.404*s+1)/(2.975*s+1));
C2 = 164.123*((5.423*s+1)/(1.153*s+1));
C3 = 350.467*((4.399*s+1)/(0.631*s+1));


L1 = C1*minreal(Gl);
L2 = C2*minreal(Gl);
L3 = C3*minreal(Gl);

figure()
nyquist(L1);
hold on
nyquist(L2);
hold on 
nyquist(L3);
xlim([-1.5 0.1]);
ylim([-1.1 1.1]);
legend("L1", "L2" , "L3" );

T1 = (L1)/(1+L1);
T2 = (L2)/(1+L2);
T3 = (L3)/(1+L3);

S1 = (1)/(1+L1);
S2 = (1)/(1+L2);
S3 = (1)/(1+L3);

MM1 = (getPeakGain(S1))^(-1);
MM2 = (getPeakGain(S2))^(-1);
MM3 = (getPeakGain(S3))^(-1);

L1 = minreal(L1);
L2 = minreal(L2);
L3 = minreal(L3);
T1_mr = minreal(T1);
T2_mr = minreal(T2);
T3_mr = minreal(T3);
 
PL1 = pole(L1);
PL2 = pole(L2);
PL3 = pole(L3);

PT1 = pole(T1_mr);
PT2 = pole(T2_mr);
PT3 = pole(T3);

%Margins (HOW TO CACULATE MM?)
[Gm1,Pm1,Wcg1,Wcp1] = margin(L1);
[Gm2,Pm2,Wcg2,Wcp2] = margin(L2);
[Gm3,Pm3,Wcg3,Wcp3] = margin(L3);



%CONST VEL
t_sim     = [50];  % Simulation time
q_init    = [10 ; pi/2];  % Initial value for q
dq_init   = [0;0];  % Initial value for dq
s_pd_of_t = [0];  % Prescribed displacement a function of time t (e.g. sin(t))
reference_signal = 'Const_vel';
reference_settings = [10 , 20 , 1 ,11];
Sim_const_vel_C1 = DOMS("Sim_Control", Components, t_sim, q_init , dq_init , s_pd_of_t, C1 , reference_signal , reference_settings);
Sim_const_vel_C2 = DOMS("Sim_Control", Components, t_sim, q_init , dq_init , s_pd_of_t, C2 , reference_signal , reference_settings);
Sim_const_vel_C3 = DOMS("Sim_Control", Components, t_sim, q_init , dq_init , s_pd_of_t, C3 , reference_signal , reference_settings);



figure()
plot(Sim_const_vel_C1.t, Sim_const_vel_C1.qt(:,1) - Components.parameters.L2.value*cos(Sim_const_vel_C1.qt(:,2)));
hold on
plot(Sim_const_vel_C2.t, Sim_const_vel_C2.qt(:,1) - Components.parameters.L2.value*cos(Sim_const_vel_C2.qt(:,2)));
hold on
plot(Sim_const_vel_C3.t, Sim_const_vel_C3.qt(:,1) - Components.parameters.L2.value*cos(Sim_const_vel_C3.qt(:,2)));
hold on
plot ([0, 1, 11, 50],[10,10,20,20])
xlim([0 50]);
ylim([9 25]);
title("Position of the load with lower $$d_{W}$$" , 'Interpreter','latex');
xlabel(str1,'Interpreter','latex');
ylabel(str4,'Interpreter','latex');
legend('Controller 1' , 'Controller 2' , 'Controller 3', 'Reference' );

%AGGR
C_aggr;
LA = C_aggr*minreal(Gl);
TA = (LA)/(1+LA);

PLA = pole(minreal(LA));
PTA = pole(minreal(TA));

figure()
nyquist(LA);
xlim([-1.5 0.5]);
ylim([-1.2 1.2]);

[Gm1,Pm1,Wcg1,Wcp1] = margin(LA);

%CONST VEL
t_sim     = [50];  % Simulation time
q_init    = [10 ; pi/2];  % Initial value for q
dq_init   = [0;0];  % Initial value for dq
s_pd_of_t = [0];  % Prescribed displacement a function of time t (e.g. sin(t))
reference_signal = 'Const_vel';
reference_settings = [10 , 20 , 1 ,11];
Sim_const_vel_C_LEAD = DOMS("Sim_Control", Components, t_sim, q_init , dq_init , s_pd_of_t, C_aggr , reference_signal , reference_settings);

figure()
plot(Sim_const_vel_C_LEAD.t, Sim_const_vel_C_LEAD.qt(:,1) - Components.parameters.L2.value*cos(Sim_const_vel_C_LEAD.qt(:,2)));
hold on
plot ([0, 1, 11, 50],[10,10,20,20])
xlim([0 50]);
ylim([9 25]);
title("Position of the load with lower $$d_{W}$$ " , 'Interpreter','latex');
xlabel(str1,'Interpreter','latex');
ylabel(str4,'Interpreter','latex');
legend("Agressive Controller" , "Reference Signal")

%%D6-------------------------------------------------------------------
 bode(Gl);

sys = Gl;

[C_pd, info_1] = pidtune(sys , "PD" ); %Add PD controller)

[C_pd, info_1] = pidtune(sys, "PD", 0.3); %Add desired crossover frequency

wc_pd = info_1.CrossoverFrequency;
PM_pd = info_1.PhaseMargin;

[C_p, info_2] = pidtune(sys, "P", 0.3); %Add desired crossover frequency

wc_p = info_2.CrossoverFrequency;
PM_p = info_2.PhaseMargin;

wc_REQ = 0.3;

K = (1)/((1+sqrt(2))*abs(evalfr(Gl, i*wc_REQ)));
tau1 = ((1+sqrt(2)))/(wc_REQ);
tau2 = (1)/((1+sqrt(2))*wc_REQ);

C_LEAD = K*(((wc_REQ)^(-1)*(1+sqrt(2))*s+1)/(((1)/(wc_REQ*(1+sqrt(2))))*s+1));

C_LEAD

L_LEAD = Gl*C_LEAD;
S_LEAD = (1)/(1+L_LEAD);

[Gm_LEAD,Pm_LEAD,Wcg_LEAD,Wcp_LEAD] = margin(L_LEAD);
MM_LEAD = (getPeakGain(S_LEAD))^(-1);

%CONST VEL
t_sim     = [50];  % Simulation time
q_init    = [10 ; pi/2];  % Initial value for q
dq_init   = [0;0];  % Initial value for dq
s_pd_of_t = [0];  % Prescribed displacement a function of time t (e.g. sin(t))
reference_signal = 'Const_vel';
reference_settings = [10 , 20 , 1 ,11];
Sim_const_vel_C_LEAD = DOMS("Sim_Control", Components, t_sim, q_init , dq_init , s_pd_of_t, C_LEAD , reference_signal , reference_settings);

figure()
plot(Sim_const_vel_C_LEAD.t, Sim_const_vel_C_LEAD.qt(:,1) - Components.parameters.L2.value*cos(Sim_const_vel_C_LEAD.qt(:,2)));
hold on
plot ([0, 1, 11, 50],[10,10,20,20])
xlim([0 50]);
ylim([9 25]);
title("Position of the load with lower $$d_{W}$$" , 'Interpreter','latex');
xlabel(str1,'Interpreter','latex');
ylabel(str4,'Interpreter','latex');
legend("Lead Controller" , "Reference Signal")

%with lower $$d_{W}$$

%%NOTCH FILTER
bode(Gl);

wc_REQ_2 = 2;

Gl_ALT = zpk(Gl);

wn_1 = sqrt(1.225);
wn_2 = (10)*wn_1;
beta_1 = 0.5*((0.5199)/(wn_1));
beta_2 = 0.9; %TUNED

N = (((s^2)/((wn_1)^2)) + ((2*beta_1)/(wn_1))*s +1)/(((s^2)/((wn_2)^2)) + ((2*beta_2)/(wn_2))*s +1);
L_N = Gl*N;


alpha = 1.84;
tau_1 = (alpha)/(wc_REQ_2);
tau_2 = (1)/(alpha*wc_REQ_2);

C_L = ((0.5*alpha^2*s+alpha)/(0.5*s + alpha));

syms K_2

fraq = Gl*C_L*N ;
eqn = K_2*abs(evalfr(fraq, i*wc_REQ_2)) == 1 ;
K2_A = solve(eqn, K_2);
K2 = 36893488147419103232/8583778407751303
C_L_K = K2 * C_L;
L_2 = Gl*C_L_K*N;

%ANALYSIS
L_1 = L_LEAD;
L_2;

T_1 = (L_1)/(1+L_1);
T_2 = (L_2)/(1+L_2);

figure()
nyquist(L_1);
hold on
nyquist(L_2);
xlim([-1.5 0.1]);
ylim([-1.1 1.1]);
legend("Lead Controller", "Lead-Notch Controller");

P_L_1 = pole(L_1);
P_L_2 = pole(L_2);
P_T_1 = pole(T_1);
P_T_2 = pole(T_2);

C_L_1 = C_LEAD;
C_L_2 = C_L_K*N ;

%CONST VEL
t_sim     = [50];  % Simulation time
q_init    = [10 ; pi/2];  % Initial value for q
dq_init   = [0;0];  % Initial value for dq
s_pd_of_t = [0];  % Prescribed displacement a function of time t (e.g. sin(t))
reference_signal = 'Const_vel';
reference_settings = [10 , 20 , 1 ,11];
Sim_const_vel_C_L_1 = DOMS("Sim_Control", Components, t_sim, q_init , dq_init , s_pd_of_t, C_L_1 , reference_signal , reference_settings);
Sim_const_vel_C_L_2 = DOMS("Sim_Control", Components, t_sim, q_init , dq_init , s_pd_of_t, C_L_2 , reference_signal , reference_settings);


figure()
plot(Sim_const_vel_C_L_1.t, Sim_const_vel_C_L_1.qt(:,1) - Components.parameters.L2.value*cos(Sim_const_vel_C_L_1.qt(:,2)));
hold on
plot(Sim_const_vel_C_L_2.t, Sim_const_vel_C_L_2.qt(:,1) - Components.parameters.L2.value*cos(Sim_const_vel_C_L_2.qt(:,2)));
hold on
plot ([0, 1, 11, 50],[10,10,20,20])
xlim([0 50]);
ylim([9 25]);
title("Position of the load with lower$$d_{W}$$ " , 'Interpreter','latex');
xlabel(str1,'Interpreter','latex');
ylabel(str4,'Interpreter','latex');
legend('Lead Controller' , 'Lead-Notch Controller' ,'Reference' );


% alpha = 2 + sqrt(3);
% 
% K2 = (1)/((alpha)*abs(evalfr(Gl, i*wc_REQ_2)));
% 
% C_LEAD_2 = K2*((0.5*alpha^2*s+alpha)/(0.5*s + alpha));


% fraq = N*C_LEAD_2*Gl;
% syms K_3 ;
% eqn = K_3*abs(evalfr(fraq, i*wc_REQ_2)) == 1 ;
% K3 = solve(eqn, K_3);
% K3_D = 36893488147419103232/6034345262679891 ; %value in decimals
%  
%  SYSTEM = K3_D*fraq ;
% K2 = (1)/((alpha)*abs(evalfr(Gl, i*wc_REQ_2)));
% 
% 
% arg_N = angle(evalfr(N, i*wc_REQ_2))*(180/(pi)); %DEGREES

% with lower $$d_{W}$$

