# Hypersonics
Purdue AAE 537 - Hypersonics Final Project Code

# Method 1
% Nozzle Contour, Rapid Method for Plug Nozzle Design
clc, clear all, close all

% Set Inputs
AR = 2;           % Expasion Area Ratio
eta_b = 0.01;     % Truncation Parameter (0=PerfectExpansion)
t_diam = 0.2;     % Throat Height (UseWhateverUnits)
step_size = 100;  % Model Fidelity
gam = 1.4;        % Ratio Specific Heats
 
% Derived Geometry
A_t = pi*((t_diam/2)^2);    % Throat Area - Assumed Circle
A_e = AR*A_t;               % Exit Area
re = sqrt(A_e/pi);          % Exit Radius

% Set Matrices
alpha = zeros(1,step_size);
A = zeros(1,step_size);
l_nd = zeros(1,step_size);
r_nd = zeros(1,step_size);
x_nd = zeros(1,step_size);
y_nd = zeros(1,step_size);
lr = zeros(1,step_size);
r = zeros(1,step_size);
xr = zeros(1,step_size);
yr = zeros(1,step_size);
mach_numbers = zeros(1,step_size);
mach_angles = zeros(1,step_size);

% Exit Mach Solver (Copied from MST Thesis!)
f1 = 2/(gam+1);
f2 = (gam-1)/2;
f3 = (gam+1)/(gam-1);
M0 = 4;
M = M0;
Function = (1/(M^2))*((f1*(1+(f2*(M^2))))^f3)-(AR^2);
Tolerance = 0.0001;
maxits = 200;
J = abs(Function);
i = 1;
while J>Tolerance && i<=maxits
    Function = (1/(M^2))*((f1*(1+(f2*(M^2))))^(f3))-(AR^2);
    dFunction = (f3/(M^2))*((f1*(1+(f2*(M^2))))^(f3-1))*2*f2*f1*M + (-2/(M^3))*((f1*(1+f2*(M^2)))^(f3));
    y = M-(Function/dFunction);
    M = y;
    Function = (1/(M^2))*((f1*(1+(f2*(M^2))))^(f3))-(AR^2);
    J = abs(Function);
    i = i+1;
end

% Use P-M to find Max Turn Angle from Exit Mach
pm1 = sqrt((gam+1)/(gam-1));
pm2 = (gam-1)/(gam+1);
pm3 = M^2-1;           
max_turn_angle = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
turn_steps2 = linspace(0,max_turn_angle,step_size);

% Contour Loop
for i = 1:step_size
    
    % Calculate Mach and PM Angles
    dummy_angle = 0;
    mach_numbers(i) = 1;
    while (dummy_angle < turn_steps2(i))
        pm1 = sqrt((gam+1)/(gam-1));
        pm2 = (gam-1)/(gam+1);
        pm3 = mach_numbers(i)^2-1;           
        dummy_angle = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
        mach_numbers(i) = mach_numbers(i) + 0.00001;
    end
    
    % Calcualte Mach and Alpha Angles
    mach_angles(i) = asind(1/mach_numbers(i));
    alpha(i) = turn_steps2(end) - turn_steps2(i) + mach_angles(i);
    
    % Calculate Streamtube Areas
    aa = (2+(gam-1)*(mach_numbers(i)^2));
    bb = (2+(gam-1)*(mach_numbers(1)^2));
    cc = (gam+1)/(2*(gam-1));
    A(i) = ( mach_numbers(1) / mach_numbers(i) ) * ( (aa / bb)^(cc) );
    
    % Calculate Non-Dimenional Parameters
    l_nd(i) = (1-sqrt(1-(A(i)*(1-(eta_b^2))...
        *mach_numbers(i)*(sind(alpha(i))/AR))))/sind(alpha(i));
    r_nd(i) = 1-(l_nd(i)*sind(alpha(i)));
    x_nd(i) = l_nd(i)*cosd(alpha(i));
    y_nd(i) = l_nd(i)*sind(alpha(i));
    
    % Calculate Coordinates
    lr(i) = l_nd(i)*re;
    r(i) = r_nd(i)*re;
    xr(i) = x_nd(i)*re*10;
    yr(i) = y_nd(i)*re*10;
    
end

figure(1)
plot(xr,yr)
xlabel('Nozzle Length')
ylabel('Nozzle Height')
grid on


# Method 2
% Nozzle Contour, Approximate Method from Characteristic Line
clc, clear all, close all

% Input Dimensions
throat_angle = 25;      % Throat Inclination
throat_height = 0.2;    % Throat Height
cowl_height = 2;        % Cowl Height to Centerline
body_width = 1;         % Body Width
step_size = 100;        % Model Fidelity
gam = 1.4;              % Ratio of Specific Heats

% Set  Matrices
x = zeros(1,step_size);
y = zeros(1,step_size);
mach_numbers = zeros(1,step_size);
mach_angles = zeros(1,step_size);
local_turn = zeros(1,step_size);
area_streamtube = zeros(1,step_size);
length_to_lip = zeros(1,step_size);

% Geometric Parameters
throat_area = throat_height*body_width;
turn_steps = linspace(0,throat_angle,step_size);

% Nozzle Contour Loop
for i = 1:step_size
    
    % Calculate Mach at all stations
    dummy_angle = 0;
    mach_numbers(i) = 1;
    while (dummy_angle < turn_steps(i))
        pm1 = sqrt((gam+1)/(gam-1));
        pm2 = (gam-1)/(gam+1);
        pm3 = mach_numbers(i)^2-1;           
        dummy_angle = pm1*atand(sqrt(pm2*pm3))-atand(sqrt(pm3));
        mach_numbers(i) = mach_numbers(i) + 0.00001;
    end
    
    % Calcualte Mach and Local Angles
    mach_angles(i) = asind(1/mach_numbers(i));
    local_turn(i) = throat_angle - turn_steps(i);
    
    % Calculate Streamtube Areas
    aa = (2+(gam-1)*(mach_numbers(i)^2));
    bb = (2+(gam-1)*(mach_numbers(1)^2));
    cc = (gam+1)/(2*(gam-1));
    area_streamtube(i) = throat_area * ( mach_numbers(1) / mach_numbers(i) ) * ( (aa / bb)^(cc) );
    
    % Calculate Contour
    length_to_lip(i) = area_streamtube(i)/sind(mach_angles(i));
    step_angle = mach_angles(i)+local_turn(i);
    x(i) = length_to_lip(i)*cosd(step_angle);
    y(i) = cowl_height - length_to_lip(i)*sind(step_angle);
    
end

figure (1)
plot(x,y)
xlabel('Nozzle Length')
ylabel('Nozzle Height')
grid on


% Modify MFP for Q addition
% 
% 
% % Exit Flow Properties
% M_throat = 1;
% mdot = 50;
% R = 287;
% T_exit = 2400;
% Tt_exit = T_exit*(1 + (gam-1)/2*M_throat^2);
% Pt_exit = 1101325;
% P_exit = Pt_exit*(1 + (gam-1)/2*M_throat^2)^(-gam/(gam-1));
% v_exit = M_throat*sqrt(gam*R*T_exit);
% Rho_exit = mdot/(v_exit*throat_area);
% Rhot_exit = Rho_exit*(Tt_exit/T_exit)^(1/(gam-1));
% 
% % Max Flow Turning, Not required
% pm_max = (pi()/2)*(sqrt((gam+1)/(gam-1))-1);
% pm_max = pm_max*180/pi();
% pm_mach_throat = sqrt((gam+1)/(gam-1))*atand(sqrt((gam-1)*(M_throat^2-1)/(gam+1)))-atand(sqrt(M_throat^2-1));
% turn_max = pm_max - pm_mach_throat;
% 
% % Set Matrices
% Static_Temperature = zeros(1,step_size);
% Stagnation_Temperature = zeros(1,step_size);
% Static_Pressure = zeros(1,step_size);
% Stagnation_Pressure = zeros(1,step_size);
% 
% % Compute flow properties along contour
% % Flow Properties Loop
% % Assume isentropic, No Q loss
% % How to model Q loss?
% for i = 1:step_size
%     
%     
%     
%     
% end







