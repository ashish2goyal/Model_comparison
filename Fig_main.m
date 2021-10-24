%% Main function simulating and plotting ODE equations
function Fig_main

% define all variables that are used by other scripts as global (for example, here we call Virology_HCV_ode.m) 
global s;
global d;
global nt;
global bt;
global T_0;
global del;
global ep;
global p;
global c;

% define parameter values
bt=7*10^-9  ; % represents viral infectivity
T_0=10^7; % represents the constant number/concentration of target cells including at time t=0
d=0.004; % natural death rate of uninfected hepatocytes
s=d*T_0; % the generation rate of uninfected hepatocytes in the absence of infection to maintain T_0 conc of target cells
c=23; % representing viral clearance


figure(1)
nt=0.998 ; % represents reduction in infectivity due to antiviral treatment
ep=0.998; % represents efficacy of the drug in inhibiting the viral production
del=0.2; % represents the death rate of infected cells 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2-dim model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since the system is assumed to be in steady state at t=0 (the time of start of antiviral treatment)
% at which we also are often aware of viral loads (V_0), we additionally have the following constraints 
V_0=10^7 ; % assumed (usually known) copies/mL level in chrnoic HCV patients
I_0=bt*V_0*T_0/del; % from dI/dt equation
p=c*V_0/I_0; % from dV/dt equation

%% Simulation for defined parameters
xinit=[I_0;V_0]; % define initial values of dependent variables in the system of ODEs in the same order as they are described 
tstop=30; % specify how long the simulation needs to be run (here we assume 30 days of treatment)
tspan=[0:0.01:tstop]; 

[t,y]=ode23s(@two_dim_model,tspan,xinit); % calling ode solver


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since the system is assumed to be in steady state at t=0 (the time of start of antiviral treatment)
% at which we also are often aware of viral loads (V_0), we additionally have the following constraints (Neumann et al., 1998) 
V_0=10^7 ; % assumed (usually known) copies/mL level in chrnoic HCV patients
T_0x= c*del / (p * bt); % from dI/dt equation
I_0=c*V_0/p; % from dV/dt equation

%% Solution
lambda1= 0.5 * ( (c+del) + sqrt((c-del)^2 + 4 * (1-ep)*(1-nt)*c*del )) 
lambda2= 0.5 * ( (c+del) - sqrt((c-del)^2 + 4 * (1-ep)*(1-nt)*c*del )) 

A1 = (ep * c - lambda2)/ (lambda1- lambda2);
A2 = (1-A1);

tspana=[0:0.01:tstop];
Va = V_0 *(A1 * exp(-lambda1 * tspana) + A2 * exp(-lambda2 * tspana));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3-dim model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since the system is assumed to be in steady state at t=0 (the time of start of antiviral treatment)
% at which we also are often aware of viral loads (V_0), we additionally have the following constraints 
V_0=10^7 ; % assumed (usually known) copies/mL level in chrnoic HCV patients
T0= s/(d + bt*V_0) ; 
I_0=bt*V_0*T0/del; % from dI/dt equation
p=c*V_0/I_0; % from dV/dt equation

%% Simulation for defined parameters
xinit=[T0; I_0;V_0]; % define initial values of dependent variables in the system of ODEs in the same order as they are described 
tstop=30; % specify how long the simulation needs to be run (here we assume 30 days of treatment)
tspan=[0:0.01:tstop]; 

[t1,y1]=ode23s(@three_dim_model,tspan,xinit); % calling ode solver

%% Plotting simulation with y-axis on log scale using semilogy

h1=semilogy(t,y(:,2),'Color','black','LineStyle','-.','LineWidth',1) % plot viral loads 
hold on
semilogy(t1,y1(:,3),'Color','green','LineStyle','--','LineWidth',1)  % plot viral loads
hold on
semilogy(tspana,Va,'Color','red','LineStyle','-','LineWidth',1)  % plot viral loads 
legend('Two-dimensional model','Three dimensional model','Analytical solution')
ylim([10^1 10^8])
xlim([0 tstop])
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('Viral loads')
xlabel('time (days)')
yticks([10^1 10^2 10^3 10^4 10^5 10^6 10^7 10^8])
yticklabels({'10^1','10^2','10^3','10^4','10^5','10^6','10^7','10^8'})
txt1 = '\eta = 0.998';
txt2 = '\epsilon = 0.998';
txt3 = '\delta = 0.2/day';
text(18,10^6.6,txt1)
text(18,10^6.0,txt2)
text(18,10^5.4,txt3)





