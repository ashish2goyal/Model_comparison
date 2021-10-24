%% the ODE equation structure

function dydt= three_dim_model(t,y) % Virology_HCV_ode will be called in an odesolver

% parameters in the model
global s;
global d;
global nt; % represents reduction in infectivity
global bt; % represents viral infectivity
global T_0; % represents the constant number/concentration of target cells
global del; % represents the death rate 
global ep; % represents efficacy of the drug in inhibiting the viral production
global p; % representing viral production
global c; % representing viral clearance

%% ODE MODEL

dydt=zeros(3,1);

dydt(1) =  s - d * y(1) - (1-nt)*bt*y(1)*y(3) ; % target cells

dydt(2) = (1-nt)*bt*y(1)*y(3)-del*y(2); % HCV infected cells (I equation; y(1))

dydt(3) = (1-ep)*p*y(2) - c*y(3); % dually-infected cells (V equation; y(2))
