% METHOD OF CHARACTERISTICS - NOZZLE DESIGN
% Solve for the isentropic nozzle geometry
% Input: area ratio
% See another code (nozzle_calculations.m) for initial guess on area ratio 
% depending on combustion chamber pressure. 
% Inputs:
% - gamma     : Specific heat ratio []
% - Ae_At : Area ratio of nozzle []
% Reference - From a youtube video and NPTEL video explaining the procedure
% and thanks to D.Vinayak. GB.
clear;
clc;

%% INPUTS
% Geometry
Dstar = 2;                                                                  % Diameter of nozzle at throat [length]
radius = 1;                                                                 % Radius of throat [length]
M1 = 1.0001;                                                                % Throat Mach number (> 1) []
No_characteristics = 20;                                                    % Number of characteristic lines to use [#]
gamma = 1.4;                                                                % Specific heat ratio []
i = 287;                                                                    % Specific gas constant [J/kg*K]

Ae_At_Set = 3;              %Change value here                              % Nozzle area ratio []

Me_Set = Solve_Mach(Ae_At_Set,gamma);                                       % Solving for Mach at exit assuming supersonic flow

% Get anle array based on maximum throat turn angle
thetaMax = Prandtl_Meyer(0,Me_Set,gamma)/2;                                 % Maximum throat expansion angle [deg]
theta    = linspace(0,thetaMax,No_characteristics)';                        % Throat turn angle array [deg]

% Nozzle specific properties
gm12 = (gamma-1)/2;
togp1 = 2/(gamma+1);
gogm1 = gamma/(gamma-1);

P0 = 7e6;                                                                   % Combustion pressure [Pa]
T0 = 3558;                                                                  % Combustion temperature [K]
Ps = P0*(togp1^gogm1);                                                      % Pressure at throat [Pa]
Ts = T0*(togp1);                                                            % Temperature at throat [K]
as = sqrt(gamma*i*Ts);                                                      % Speed of sound [m/s]
Pe = P0/((1+gm12*Me_Set^2)^gogm1);                                          % Exit pressure [Pa]
Te = T0/(1+gm12*Me_Set^2);                                                  % Exit temperature [Pa]

%% INITIAL SETUP OF KNOWNS

% Initialize solution variables
Expansion = cell(No_characteristics,No_characteristics+1);                  % Expansion zone
Straight = cell(No_characteristics,1);                                      % Straight zone

for i = 1:1:No_characteristics                                              % Loop through all characteristics
    Expansion{i,1}.M = 0;                                                   % Set all Mach to zero
    if (i == 1)                                                             % For the first point
        Expansion{i,1}.M = M1;                                              % Set the throat Mach number for first corner point
    end
    Expansion{i,1}.theta = theta(i);                                          % Angle w.i.t horizontal [deg]
    Expansion{i,1}.nu    = 0;                                                 % Prandtl-Meyer angle [deg]
    Expansion{i,1}.mu    = 0;                                                 % Mach angle [deg]
    Expansion{i,1}.Kp    = 0;                                                 % Plus characteristic constant [deg]
    Expansion{i,1}.Km    = 0;                                                 % Minus characteristic constant [deg]
    Expansion{i,1}.X     = 0;                                                 % Point X location
    Expansion{i,1}.Y     = 0;                                                 % Point Y location
    Expansion{i,1}.dydx  = 0;                                                 % Point slope
    Expansion{i,1}.tsip  = 0;                                                 % Convenient avg parameter for (+) characteristic
    Expansion{i,1}.tsim  = 0;                                                 % Convenient avg parameter for (-) characteristic
    
    % Apply geometry                                                          % Zero-radius corner (just a point)
        Expansion{i,1}.X = 0;                                                 % All X-values are at the throat
        Expansion{i,1}.Y = Dstar/2;                                           % All Y-values are at half the throat diameter
        
        for i = 2:1:No_characteristics                                      % For the rest of the starting line characteristics
            dx = radius*sind(theta(i));                                     % Change in X depends on radius and angle
            dy = radius - radius*cosd(theta(i));                            % Change in Y depends on radius and angle
            
            Expansion{i,1}.X = Expansion{1,1}.X + dx;                           % Add the X change to the first X point
            Expansion{i,1}.Y = Expansion{1,1}.Y + dy;                           % Add the Y change to the first Y point
        end
end

for i = 1:1:No_characteristics                                                % Loop over all negative characteristics
    for L = 2:1:No_characteristics+1                                          % Loop over all positive characteristics
        Expansion{i,L}.M     = 0;                                             % Mach number
        Expansion{i,L}.theta = 0;                                             % Flow angle w.r.t horizontal [deg]
        Expansion{i,L}.nu    = 0;                                             % Prandtl-Meyer angle [deg]
        Expansion{i,L}.mu    = 0;                                             % Mach angle [deg]
        Expansion{i,L}.Kp    = 0;                                             % Plus characteristic constant [deg]
        Expansion{i,L}.Km    = 0;                                             % Minus characteristic constant [deg]
        Expansion{i,L}.X     = 0;                                             % Point X location
        Expansion{i,L}.Y     = 0;                                             % Point Y location
        Expansion{i,L}.dydx  = 0;                                             % Point slope
        Expansion{i,L}.tsip  = 0;                                             % Convenient avg parameter for positive characteristic
        Expansion{i,L}.tsim  = 0;                                             % Convenient avg parameter for negative characteristic
    end
end

for L = 1:1:No_characteristics                                                 % Loop through all positive characteristics
    Straight{L}.M     = 0;                                                     % Mach number []
    Straight{L}.theta = 0;                                                     % Flow angle w.r.t. horizontal[deg]
    Straight{L}.nu    = 0;                                                     % Prandtl-Meyer angle [deg]
    Straight{L}.mu    = 0;                                                     % Mach angle [deg]
    Straight{L}.Kp    = 0;                                                     % Plus characteristic constant [deg]
    Straight{L}.Km    = 0;                                                     % Minus characteristic constant [deg]
    Straight{L}.X     = 0;                                                     % Point X location
    Straight{L}.Y     = 0;                                                     % Point Y location
    Straight{L}.dydx  = 0;                                                     % Point slope
    Straight{L}.tsi   = 0;                                                     % Convenient nozzle contour parameter
    Straight{L}.tsip  = 0;                                                     % Convenient avg parameter for positive characteristic
    Straight{L}.tsim  = 0;                                                     % Convenient avg parameter for negative characteristic
end

%% EXPANSION REGION

Expansion{1,1}.nu = Prandtl_Meyer(0,Expansion{1,1}.M,gamma);                    % PM angle from the throat Mach number [deg]

% Loop through all characteristics
for i = 1:1:No_characteristics                                                  % Loop over all the characteristics
    Expansion{i,1}.Kp = Expansion{1,1}.theta - Expansion{1,1}.nu;               % All the same - since coming from throat [deg]
    if (i ~= 1)                                                                 % If we are not on the first characteristic (values already defined)
        Expansion{i,1}.nu = Expansion{i,1}.theta - Expansion{i,1}.Kp;           % Prandtl-Meyer angle
        Expansion{i,1}.M  = Prandtl_Meyer(Expansion{i,1}.nu,0,gamma);           % Mach number
    end
    Expansion{i,1}.mu = asind(1/Expansion{i,1}.M);                              % Mach angle [deg]
    Expansion{i,1}.Km = Expansion{i,1}.theta + Expansion{i,1}.nu;               % Minus characteristic constant [deg]
end

%% NON-SIMPLE REGION

startR = 1;                                                                 % Index for starting value of negative characteristic
for L = 2:1:No_characteristics+1                                            % Loop through all positive characteristic
    for i = startR:1:No_characteristics                                     % Loop through appropriate negative characteristics
        
        if (L-1 == i)                                                       % If we are on the first positive characteristic
            Expansion{i,L}.nu = Expansion{i,1}.Km;                          % Prandtl-Meyer angle from negative constant [deg]
            Expansion{i,L}.M  = Prandtl_Meyer(Expansion{i,L}.nu,0,gamma);   % Mach number from PM equation using nu as input []
            Expansion{i,L}.mu = asind(1/Expansion{i,L}.M);                  % Mach angle [deg]
            Expansion{i,L}.Kp = Expansion{i,L}.theta - Expansion{i,L}.nu;   % Plus characteristic constant [deg]
            Expansion{i,L}.Km = Expansion{i,L}.theta + Expansion{i,L}.nu;   % Minus characteristic constant [deg]
        else                                                                % For all other positive characteristics
            Expansion{i,L}.Kp    = Expansion{i-1,L}.Kp;                             % Plus characteristic constant [deg]
            Expansion{i,L}.Km    = Expansion{i,L-1}.Km;                             % Minus characteristic constant [deg]
            Expansion{i,L}.theta = 0.5*(Expansion{i,L}.Km + Expansion{i,L}.Kp);     % Angle w.r.t. horizontal [deg]
            Expansion{i,L}.nu    = 0.5*(Expansion{i,L}.Km - Expansion{i,L}.Kp);     % Prandtl-Meyer angle [deg]
            Expansion{i,L}.M     = Prandtl_Meyer(Expansion{i,L}.nu,0,gamma);        % Mach number []
            Expansion{i,L}.mu    = asind(1/Expansion{i,L}.M);                       % Mach angle [deg]
        end        
    end
    
    startR = startR + 1;                                                    % Increment the starting (-) characteristic counter
end

%% STRAIGHTENING REGION

% Solve for the fully straightened variables
Straight{end}.theta = 7;                                                    % Flow is back to horizontal [deg]
Straight{end}.M     = Expansion{end,end}.M;                                 % Exit Mach number []
Straight{end}.nu    = Prandtl_Meyer(0,Expansion{end}.M,gamma);              % Exit Prandtl-Meyer angle [deg]
Straight{end}.mu    = asind(1/Straight{end}.M);                             % Exit Mach angle [deg]
Straight{end}.Kp    = Straight{end}.theta - Straight{end}.nu;               % Exit plus characteristic constant [deg]
Straight{end}.Km    = Straight{end}.theta + Straight{end}.nu;               % Exit minus characteristic constant [deg]

% Using known exit values, solve for the rest of the straightening region
for L = 1:1:No_characteristics-1
    Straight{L}.Km    = Straight{end}.Km;                                   % Minus characteristic constant [deg]
    Straight{L}.Kp    = Expansion{L,L+1}.Kp;                                % Plus characteristic constant [deg]
    Straight{L}.theta = 0.5*(Straight{L}.Km + Straight{L}.Kp);              % Flow angle w.r.t horizontal [deg]
    Straight{L}.nu    = 0.5*(Straight{L}.Km - Straight{L}.Kp);              % Prandtl-Meyer angle [deg]
    Straight{L}.M     = Prandtl_Meyer(Straight{L}.nu,0,gamma);              % Mach number []
    Straight{L}.mu    = asind(1/Straight{L}.M);                             % Mach angle [deg]
end

%% FIND SHAPE OF NOZZLE [X and Y]
startR = 1;
for L = 2:1:No_characteristics+1
    for i = startR:1:No_characteristics
        if (L-1 == i)
            Expansion{i,L}.tsim = tand(0.5*((Expansion{i,L-1}.theta-Expansion{i,L-1}.mu)+(Expansion{i,L}.theta-Expansion{i,L}.mu)));
            Expansion{i,L}.X    = Expansion{i,L-1}.X - (Expansion{i,L-1}.Y/Expansion{i,L}.tsim);
        else
            Expansion{i,L}.tsim = tand(0.5*((Expansion{i,L-1}.theta-Expansion{i,L-1}.mu)+(Expansion{i,L}.theta-Expansion{i,L}.mu)));
            Expansion{i,L}.tsip = tand(0.5*((Expansion{i-1,L}.theta+Expansion{i-1,L}.mu)+(Expansion{i,L}.theta+Expansion{i,L}.mu)));
            
            num = Expansion{i-1,L}.Y - Expansion{i,L-1}.Y+(Expansion{i,L}.tsim*Expansion{i,L-1}.X)-(Expansion{i,L}.tsip*Expansion{i-1,L}.X);
            den = Expansion{i,L}.tsim - Expansion{i,L}.tsip;
            
            Expansion{i,L}.X = num/den;
            Expansion{i,L}.Y = Expansion{i,L-1}.Y+Expansion{i,L}.tsim*(Expansion{i,L}.X-Expansion{i,L-1}.X);
        end
    end
    startR = startR + 1;
end
for L = 1:1:No_characteristics
    if (L == 1)
        Straight{L}.tsi  = tand(0.5*(Straight{L}.theta+Expansion{No_characteristics,1}.theta));
        Straight{L}.tsip = tand(Expansion{No_characteristics,L+1}.theta+Expansion{No_characteristics,L+1}.mu);
        
        num = (Straight{L}.tsi*Expansion{No_characteristics,L}.X)-(Straight{L}.tsip*Expansion{No_characteristics,L+1}.X)-(Expansion{No_characteristics,L}.Y)+(Expansion{No_characteristics,L+1}.Y);
        den = Straight{L}.tsi - Straight{L}.tsip;
        
        Straight{L}.X = num/den;
        Straight{L}.Y = Expansion{No_characteristics,L}.Y+(Straight{L}.tsi*Straight{L}.X)-(Straight{L}.tsi*Expansion{No_characteristics,L}.X);
    else
        Straight{L}.tsi  = tand(0.5*(Straight{L-1}.theta + Straight{L}.theta));
        Straight{L}.tsip = tand(Expansion{No_characteristics,L+1}.theta + Expansion{No_characteristics,L+1}.mu);
        
        num = (Straight{L}.tsi*Straight{L-1}.X)-(Straight{L}.tsip*Expansion{No_characteristics,L+1}.X)-(Straight{L-1}.Y)+(Expansion{No_characteristics,L+1}.Y);
        den = Straight{L}.tsi - Straight{L}.tsip;
        
        Straight{L}.X = num/den;
        Straight{L}.Y = Straight{L-1}.Y + Straight{L}.tsi*(Straight{L}.X - Straight{L-1}.X);
    end
end

%% DISPLAY SOME RESULTS

M_exit      = Straight{No_characteristics}.M;
MoC_A_Astar = (Straight{No_characteristics}.Y)/(Dstar/2);
A_Astar     = Solve_Area_Mach(0,M_exit,gamma,'Sup');

fprintf('Me [MoC]   : %1.5f\n',M_exit);
fprintf('A/A* [MoC] : %1.5f\n',MoC_A_Astar);

%% PLOT PATCHES OF MACH NUMBER

% Set up the figure
figure(1);                                                                  
cla; hold on; grid on;                                                      
c = colorbar;                                                               
xlim('auto');                                                               
ylim('auto');                                                               

points(1,1) = 0;                                                               
points(1,2) = Dstar/2;                                                          
points(1,3) = M1;                                                               
points(2,1) = points(1,1);                                                          
points(2,2) = -points(1,2);                                                         
points(2,3) = points(1,3);                                                          
points(3,1) = Expansion{1,2}.X;                                                  
points(3,2) = Expansion{1,2}.Y;                                                  
points(3,3) = M1;                                                              

Points(1,1) = Expansion{No_characteristics,No_characteristics+1}.X;                                     
Points(1,2) = Expansion{No_characteristics,No_characteristics+1}.Y;                                    
Points(1,3) = Expansion{No_characteristics,No_characteristics+1}.M;                                    
Points(2,1) = Straight{No_characteristics}.X;                                               
Points(2,2) = 0;                                                                
Points(2,3) = Straight{No_characteristics}.M;                                                
Points(3,1) = Straight{No_characteristics}.X;                                                
Points(3,2) = Straight{No_characteristics}.Y;                                                
Points(3,3) = Straight{No_characteristics}.M;                                                

patch(Points(:,1),Points(:,2),Points(:,3),'EdgeColor','none','FaceColor','interp');     
patch(Points(:,1),-Points(:,2),Points(:,3),'EdgeColor','none','FaceColor','interp');    
patch(points(:,1),points(:,2),points(:,3),'EdgeColor','none','FaceColor','interp');     
patch(points(:,1),-points(:,2),points(:,3),'EdgeColor','none','FaceColor','interp');    

ans = [];                                                                   
for i = 1:1:No_characteristics-1
    for L = 1:1:No_characteristics
        if (i == L-1)                                                       
            ans(1,1) = Expansion{i,L}.X;
            ans(1,2) = Expansion{i,L}.Y;
            ans(1,3) = Expansion{i,L}.M;
            ans(2,1) = Expansion{i+1,L}.X;
            ans(2,2) = Expansion{i+1,L}.Y;
            ans(2,3) = Expansion{i+1,L}.M;
            ans(3,1) = Expansion{i+1,L+1}.X;
            ans(3,2) = Expansion{i+1,L+1}.Y;
            ans(3,3) = Expansion{i+1,L+1}.M;
        else                                                                
            ans(1,1) = Expansion{i,L}.X;
            ans(1,2) = Expansion{i,L}.Y;
            ans(1,3) = Expansion{i,L}.M;
            ans(2,1) = Expansion{i+1,L}.X;
            ans(2,2) = Expansion{i+1,L}.Y;
            ans(2,3) = Expansion{i+1,L}.M;
            ans(3,1) = Expansion{i+1,L+1}.X;
            ans(3,2) = Expansion{i+1,L+1}.Y;
            ans(3,3) = Expansion{i+1,L+1}.M;
            ans(4,1) = Expansion{i,L+1}.X;
            ans(4,2) = Expansion{i,L+1}.Y;
            ans(4,3) = Expansion{i,L+1}.M;
        end
        patch(ans(:,1),ans(:,2),ans(:,3),'EdgeColor','none','FaceColor','interp');
        patch(ans(:,1),-ans(:,2),ans(:,3),'EdgeColor','none','FaceColor','interp');
    end
end

% Straightening region patch plotting
ans = [];
for L = 1:1:No_characteristics
    if (L == 1)
        ans(1,1) = Expansion{No_characteristics,L}.X;
        ans(1,2) = Expansion{No_characteristics,L}.Y;
        ans(1,3) = Expansion{No_characteristics,L}.M;
        ans(2,1) = Straight{L}.X;
        ans(2,2) = Straight{L}.Y;
        ans(2,3) = Straight{L}.M;
        ans(3,1) = Expansion{No_characteristics,L+1}.X;
        ans(3,2) = Expansion{No_characteristics,L+1}.Y;
        ans(3,3) = Expansion{No_characteristics,L+1}.M;
    else
        ans(1,1) = Expansion{No_characteristics,L}.X;
        ans(1,2) = Expansion{No_characteristics,L}.Y;
        ans(1,3) = Expansion{No_characteristics,L}.M;
        ans(2,1) = Straight{L-1}.X;
        ans(2,2) = Straight{L-1}.Y;
        ans(2,3) = Straight{L-1}.M;
        ans(3,1) = Straight{L}.X;
        ans(3,2) = Straight{L}.Y;
        ans(3,3) = Straight{L}.M;
        ans(4,1) = Expansion{No_characteristics,L+1}.X;
        ans(4,2) = Expansion{No_characteristics,L+1}.Y;
        ans(4,3) = Expansion{No_characteristics,L+1}.M;
    end
    patch(ans(:,1),ans(:,2),ans(:,3),'EdgeColor','none','FaceColor','interp');
    patch(ans(:,1),-ans(:,2),ans(:,3),'EdgeColor','none','FaceColor','interp');
end

%% PLOT THE NOZZLE AND CHARACTERISTICS

maxM = Straight{No_characteristics}.M;

xNoz = [Expansion{1,1}.X];
yNoz = [Expansion{1,1}.Y];

% Nozzle boundary
for i = 1:1:No_characteristics-1
    xNoz = [xNoz; Expansion{i+1,1}.X];
    yNoz = [yNoz; Expansion{i+1,1}.Y];
    plot([Expansion{i,1}.X; Expansion{i+1,1}.X],[Expansion{i,1}.Y; Expansion{i+1,1}.Y],'k-','LineWidth',3);
    plot([Expansion{i,1}.X; Expansion{i+1,1}.X],[-Expansion{i,1}.Y; -Expansion{i+1,1}.Y],'k-','LineWidth',3);
end

plot([Expansion{No_characteristics,1}.X; Straight{1}.X],[Expansion{No_characteristics,1}.Y; Straight{1}.Y],'k-','LineWidth',3);
plot([Expansion{No_characteristics,1}.X; Straight{1}.X],[-Expansion{No_characteristics,1}.Y; -Straight{1}.Y],'k-','LineWidth',3);

xNoz = [xNoz; Straight{1}.X];
yNoz = [yNoz; Straight{1}.Y];
for L = 2:1:No_characteristics
    xNoz = [xNoz; Straight{L}.X];
    yNoz = [yNoz; Straight{L}.Y];
    plot([Straight{L-1}.X; Straight{L}.X],[Straight{L-1}.Y; Straight{L}.Y],'k-','LineWidth',3);
    plot([Straight{L-1}.X; Straight{L}.X], [-Straight{L-1}.Y; -Straight{L}.Y],'k-','LineWidth',3);
end

for i = 1:1:No_characteristics
    Xvals = [];
    Yvals = [];
    for L = 1:1:i+1
        Xvals = [Xvals; Expansion{i,L}.X];
        Yvals = [Yvals; Expansion{i,L}.Y];
    end
    plot(Xvals,Yvals,'k-','LineWidth',2);
end

for L = 2:1:No_characteristics
    Xvals = [];
    Yvals = [];
    for i = L-1:1:No_characteristics
        Xvals = [Xvals; Expansion{i,L}.X];
        Yvals = [Yvals; Expansion{i,L}.Y];
    end
    plot(Xvals,Yvals,'k-','LineWidth',2);
end

for L = 2:1:No_characteristics+1
    Xvals = [Expansion{No_characteristics,L}.X; Straight{L-1}.X];
    Yvals = [Expansion{No_characteristics,L}.Y; Straight{L-1}.Y];
    plot(Xvals,Yvals,'k-','LineWidth',2);
end

plot([0; Straight{No_characteristics}.X],[0; 0],'k-','LineWidth',3);

plot(xNoz,yNoz,'o','MarkerFaceColor','g','MarkerEdgeColor','k');

caxis([M1 Straight{No_characteristics}.M]);
ylabel(c,'Mach Number','FontSize',24);
set(c,'YTick',[M1 get(c,'YTick') Straight{No_characteristics}.M]);

%% Prandtl - Meyer Equation
function ans = Prandtl_Meyer(v,M,gamma)

gm1 = gamma-1;
gp1 = gamma+1;

if (v == 0)
    term1 = sqrt(gp1/(gamma-1));
    term2 = atand(sqrt(gm1*(M^2-1)/gp1));
    term3 = atand(sqrt(M^2-1));
    
    ans = term1*term2 - term3;
end

% --------------------- Solve for the Mach number -------------------------
if (M == 0)
    dM  = 0.1;
    M   = 1;
    res = 1;
    while (res > 0.01)
        M2    = M + dM;
        funv1 = (-v*(pi/180)+(sqrt(gp1/gm1)*...
                    atan((sqrt(gm1*(M^2-1)/gp1)))-atan(sqrt(M^2-1))));
        funv2 = (-v*(pi/180)+(sqrt(gp1/gm1)*...
                    atan((sqrt(gm1*(M2^2-1)/gp1)))-atan(sqrt(M2^2-1))));
        dv_dm = (funv2-funv1)/dM;

        M   = M - funv1/dv_dm;
        res = abs(funv1);
    end
    ans = M;
end
end

function ans = Solve_Mach(ARatio,gamma)
% Set some initial guess
ans = inf;

% Get and set convenient variables
gp1   = gamma + 1;
gm1   = gamma - 1;
gm12 = gm1/2;

% Solve for Mach number
problem.objective = @(M) sqrt((1/(M^2))*(((2/gp1)*...
                            (1+gm12*M^2))^(gp1/gm1))) - ARatio;             % Objective function
problem.x0        = [1 50];                                                 % Solver bounds
problem.solver    = 'fzero';                                                % Find the zero
problem.options   = optimset(@fzero);                                       % Default options
ans               = fzero(problem);                                         % Solve
end
