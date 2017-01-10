% SiStER_Input_File


% DURATION OF SIMULATION AND FREQUENCY OF OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%
Nt=200; % max number of time iterations
dt_out=10; % output files every "dt_out" iterations

% DOMAIN SIZE AND GRIDDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsize=800e3;
ysize=150e3;
% gridding- from 0 to GRID.x(1), grid size is GRID.dx(1)
% from GRID.x(1) to GRID.x(2), grid size is GRID.dx(1) etc...
% same for y
GRID.dx(1)=2000;
GRID.x(1)=380e3;
GRID.dx(2)=1000;
GRID.x(2)=420e3;
GRID.dx(3)=1600;
GRID.x(3)=800e3;
GRID.dy(1)=1600;
GRID.y(1)=46e3;
GRID.dy(2)=1000;
GRID.y(2)=56e3;
GRID.dy(3)=1600;
GRID.y(3)=150e3;


% LAGRANGIAN MARKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mquad=6; % number of markers in the smallest quadrant
Mquad_crit=4; % minimum number of markers allowed in smallest quadrant (for reseeding)

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nphase=4; % number of phases

% phase 1
GEOM(1).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(1).top=0;
GEOM(1).bot=50e3;

% phase 2
GEOM(2).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(2).top=50e3;
GEOM(2).bot=100e3;

% phase 3
GEOM(3).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(3).top=100e3;
GEOM(3).bot=150e3;

% phase 4http://cdn.twentytwowords.com/wp-content/uploads/9hus.jpg?b24a88
GEOM(4).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(4).top=150e3;
GEOM(4).bot=150e3;

% FZ GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% triggers temp structure and bathymetric step for n oceanic lithosphere
% blocks
GEOM(2).fzSW=1;                                 % fracture zone switch
GEOM(2).fz=[0.45 0.475 0.5 0.525 0.55 1.0];           % fraction of domain for each lithosphere block (n entries)
GEOM(2).age=3.154e13*[5 4 3 2 1 0.1];%[2.5 2.4 2.3 2.2 2.1 0.1];  % age of each block [Ma]*convert
GEOM(2).step=2*3.1e-5*3300*1280/2300*sqrt(3*abs(GEOM(2).age-max(GEOM(2).age))/(3300*1000)/pi);
%disp(num2str(GEOM(2).step(1:2)))
%340*sqrt(abs(GEOM(2).age-max(GEOM(2).age))/3.154e13);               % bathymetric step size
GEOM(2).FZseed=1e3*[5 100];   %weak seed dimensions [width depth] km

GEOM(2).fzstrong=[0 0 0 0 1];

% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creep laws of the form: pre^(-1/n)*epsII^((1-n)/n)*exp(E/(nRT))
% harmonically averaging diffusion creep, dislocation creep 
% (and plastic creep to simulate brittle failure)

% phase 1
MAT(1).phase=1;
% density parameters
MAT(1).rho0=1.025e3;     
MAT(1).alpha=0;
% elasticity 
MAT(1).G=1e18;
% diffusion creep parameters
MAT(1).pre_diff=.5/1e18;
MAT(1).Ediff=0;
MAT(1).ndiff=1;
% dislocation creep parameters
MAT(1).pre_disc=.5/1e18;
MAT(1).Edisc=0;
MAT(1).ndisc=1;
% plasticity
MAT(1).mu=0.6;
MAT(1).Cmax=1e6;
MAT(1).Cmin=0.01e6;
MAT(1).ecrit=0.1;


% phase 2
MAT(2).phase=2;
% density parameters
%MAT(2).rho0=1e3;
%MAT(2).alpha=0;
MAT(2).rho0=3330;
MAT(2).alpha=3.1e-5;
% elasticity 
MAT(2).G=30e9;
% diffusion creep parameters
%MAT(2).pre_diff=.5/1e18;
%MAT(2).Ediff=0;
%MAT(2).ndiff=1;
% dislocation creep parameters
%MAT(2).pre_disc=.5/1e18;
%MAT(2).Edisc=0;
%MAT(2).ndisc=1;
% diffusion creep parameters
MAT(2).pre_diff=0.5/1e40;
MAT(2).Ediff=0;
MAT(2).ndiff=1;
% dislocation creep parameters
MAT(2).pre_disc=1e-3;
MAT(2).Edisc=540000;
MAT(2).ndisc=1;
% plasticity
MAT(2).mu=0.6;
MAT(2).Cmax=1e6;
MAT(2).Cmin=0.1e6;
MAT(2).ecrit=0.5;


% phase 3
MAT(3).phase=3;
% density parameters
MAT(3).rho0=3330;
MAT(3).alpha=3.1e-5;
% elasticity 
MAT(3).G=30e9;
% diffusion creep parameters
MAT(3).pre_diff=0.5/1e40;
MAT(3).Ediff=0;
MAT(3).ndiff=1;
% dislocation creep parameters
MAT(3).pre_disc=1e-3;
MAT(3).Edisc=540000;
MAT(3).ndisc=1;
% plasticity
MAT(3).mu=0.6;
MAT(3).Cmax=1e6;
MAT(3).Cmin=0.1e6;
MAT(3).ecrit=0.5;

% phase 4
MAT(4).phase=4;
% density parameters
MAT(4).rho0=3330;
MAT(4).alpha=3.1e-5;
% elasticity 
MAT(4).G=30e9;
% diffusion creep parameters
MAT(4).pre_diff=0.5/1e20;
MAT(4).Ediff=0;
MAT(4).ndiff=1;
% dislocation creep parameters
MAT(4).pre_disc=0.5/1e20;
MAT(4).Edisc=0;
MAT(4).ndisc=1;
% plasticity
MAT(4).mu=0.6;
MAT(4).Cmax=1e6;
MAT(4).Cmin=0.1e6;
MAT(4).ecrit=0.5;

% ADDITIONAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.YNElast=1; % elasticity on (1) or off (0)
PARAMS.YNPlas=1; % plasticity on (1) or off (0)
PARAMS.epsII_from_stress=0; % get strain rate from stresses (1, default) or from velocity field (0)
PARAMS.tau_heal=1e40; % healing time for plasticity (s)
PARAMS.gx=0; % gravity along x
PARAMS.gy=9.8; % gravity along y
PARAMS.fracCFL=0.5; % distance by which a marker is allowed to move over a time step, expressed as a fraction of the smallest cell size
PARAMS.R=8.314; % gas constant
PARAMS.etamax=1e28; % maximum viscosity
PARAMS.etamin=1e16; % minimum viscosity
PARAMS.Tsolve=1; % yes (1) or no (0) solve for temperature
% initial temperature profile, polynomial with depth 
% T = a0 + a1*y+a2*y^2+a3*y^3+amp*sin(2*pi*X/lam)
% (make sure it matches the BCs)
PARAMS.a0=0;
PARAMS.a1=0;
PARAMS.a2=0;
PARAMS.a3=1000/(30e3)^3;
PARAMS.amp=0; % amplitude of sinusoidal perturbation
PARAMS.lam=1; % wavelength of sinusoidal perturbation
PARAMS.ynTreset=1; % if ==1, reset T=T0 where im==1 (sticky layer)
PARAMS.T0=20;
% reference values for the constant diffusivity thermal solver
% (kappa = kref / (rhoref*cpref))
PARAMS.rhoref=MAT(2).rho0; 
PARAMS.kref=2.66;
PARAMS.cpref=1000;


% TOPOGRAPHY EVOLUTION (interface between rock and sticky air/water layer)
PARAMS.Ntopo_markers=1000; % number of markers in marker chain tracking topography
PARAMS.YNSurfaceProcesses=0; % surface processes (diffusion of topography) on or off
PARAMS.topo_kappa=1e-9; % diffusivity of topography (m^2/s)

% Solver iterations
PARAMS.Npicard_min=1; % minimum number of Picard iterations per time step
PARAMS.Npicard_max=50; % maximum number of Picard iterations per time step
PARAMS.conv_crit_ResL2=1e-10;
PARAMS.pitswitch=30; % number of Picard iterations at which the solver switches to quasi-Newton
% convergence is assumed if a fraction of the domain smaller 
% than PARAMS.conv_crit2 * domain size has more than relative
% change of PARAMS.conv_crit1 in strain rate and velocity field


% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pressure
PARAMS.p0cell=0; % pressure in the top-left corner of the domain (anchor point)


% flow

% boundary conditions
% entries in BC correspond to
% 1/ rollers? (1=yes, 0=no)
% 2/ type of velocity normal to boundary (0=constant, 1=triangular profile, set up for top and bottom only)
% 3/ value of normal velocity (or max. value along the profile, if 2/ is 1)
Uin=0;%-6.634e-6;

BC.top=[0 0 0];%1.0563e-10];
BC.bot=[1 3 0];%4e-5]; %-1.0563e-10];
BC.left=[1 2 0]; %-3.1688e-10];
BC.right=[1 2 0]; %3.1688e-10];
PARAMS.BalanceStickyLayer=0; % if set to 1, the code will reset the inflow 
% / outflow BCs to balance the inflow / outflow of sticky layer material,
% and rock separately, based on the position of the sticky layer / air
% interface


% Mixed BC; Only matters if BC(1,2)=2
BC.left_w  = [1 0 0]; % Water
BC.left_l  = [1 0 0.5*Uin]; % Lithospheere
BC.left_m  = [1 0 0]; % Mantle

BC.right_w = [1 0 0]; % Water
BC.right_l = [1 0 -0.5*Uin]; % Lithosphere
BC.right_m = [1 0 0]; % Mantle

%BC.bot_profile= FIX THIS LATER

BC.sliding_lith  = [1 3]; % [a b]; if a=1, applies lith. BC to b elements
                          % above and below lithosphere. This is to allow
                          % lithosphere sliders to move without sliders in
                          % the water and mantle 
BC.mech_lith = 0;         % Use viscosity to define lithosphere vs layer 2 
BC.min_lith_visc = 1e20; % Viscosities above this are considered to be
                           % lithosphere for mixed BC purposes
                           
% BC.sticky_bot_left = 20e3;
% BC.sticky_bot_right = 20e3;
% BC.lith_bot_left = 250e3;
% BC.lith_bot_right = 250e3;
% 
% BC.p_lith_bot=8.200e9;
% BC.p_lith_top=0;

BC.INJ.on=0;

% thermal 

% entries in BCtherm correspond to
% 1/ type? (1=Dirichlet, 0=Neumann)
% 2/ value
BCtherm.top=[1 0];
BCtherm.bot=[1 1365];
BCtherm.left=[2 0];
BCtherm.right=[2 0];
