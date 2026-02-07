% Eclipsing Binary Stars Light Curve Simulation

%{
Major Assumptions:
- Uniform Luminosity Distribution (No limb darkening, luminosity/projected area = const)
- Circular Orbit (e = 0)
- Celestial bodies m1 and m2 is perfectly spherical with projected areas being perfect circles
%}


% General Scientific Constants (DO NOT CHANGE)-
AU = 1.5e11; % Astronomical Unit
G = 6.67e-11; % Gravitational Constant
% --------------------------------------------

% Astrophysics (DO NOT CHANGE)----------------
M_Sol = 1.99e30; % Solar Mass
L_Sol = 3.9e26; % Solar Luminosity
R_Sol = 6.96e8; % Solar Radii
yr = 365 * 24 * 60 * 60; % Years in seconds
% --------------------------------------------


% Code ---------------------------------------

% INPUT PARAMETERS 
m1 = .5; % Solar Masses
m2 = .25; % Solar Masses
r1 = .6; % Solar Radii
r2 = .5; % r2 < r1 Solar Radii
L1 = 1; % Solar Luminosity
L2 = 0.5; % Solar Luminosity
P = 0.154; % Orbital Period in Days
i = 85; % Inclination Angle in Deg

% Calculations (DO NOT CHANGE)
A1 = pi * r1^2; % Area of Star 1 in R_Solar^2
A2 = pi * r2^2; % Area of Star 2 in R_Solar^2
A_total = A1 + A2; % Total Area R_Solar^2
P = P * 24 * 60 * 60;
sma = ((P/yr)^2 * (m1 + m2))^(1/3) * AU; % SMA in m
v_rel = sqrt((G*(m1 + m2) * M_Sol)/sma); % Relative Orbital Velocity in m/s
b_h = sma*cos(i * pi/180)/R_Sol; % Impact Parameter Physical Length in Solar Radii
w = (2 * pi)/P; % Angular Velocity of system in rad/s
phi = asin(sqrt((r1 + r2)^2 - b_h^2) * R_Sol/sma); % Phase Shift for graph in rad
t_total = P/pi * asin(sqrt((r1 + r2)^2 - b_h^2) * R_Sol/sma); % Transit Duration in s
A1 = pi * r1^2; % Projected Area of Star 1 in R_Sol^2
A2 = pi * r2^2; % Projected Area of Star 2 in R_Sol^2
L_total = L1 + L2; % Total Luminosity of System in L_Sol


% Boundary Limits (DO NOT CHANGE)
P_min = sqrt((((r1+r2)*R_Sol)/AU)^3/(m1+m2)) * yr; % Minimum Orbital Period in Seconds
i_min = acos(((r1 + r2)*R_Sol)/sma) * (180/pi); % Minimum Inclination for eclipse in Deg
i_grazing = acos(((r1 - r2)*R_Sol)/sma) * (180/pi); % Maximum Inclination for grazing eclipse

% Equations (DO NOT CHANGE)
v = @(t) v_rel*cos(w*t - phi); % Radial Velocity as a function of t
p = @(t) ((v_rel/w)/R_Sol) * sin(w*t - phi); % Starting Position of center of Star 2, = 0 at first contact
d = @(t) sqrt(b_h^2 + p(t)^2);


% Projected Area Eclipsed in R_Sol^2 as a function of time (DO NOT CHANGE)
A_c = @(t) r1^2 * acos((d(t)^2 + r1^2 - r2^2)/(2*d(t)*r1)) + r2^2*acos((d(t)^2 + r2^2 - r1^2)/(2*d(t)*r2)) - 0.5 * sqrt((d(t)^2 -(r2 - r1)^2) * ((r1 + r2)^2 - d(t)^2));


% Times (DO NOT CHANGE)
t1 = 0; % First Contact


% Primary Eclipse Graphs (DO NOT CHANGE)
L_PE1 = @(t) (L2 + ((A1 - A_c(t))/A1) * L1);
L_PE2 = @(t) (L2 + ((A1 - A2)/A1) * L1);


% Secondary Eclipse Graphs (DO NOT CHANGE)
L_SE1 = @(t) (L1 + ((A2 - A_c(t))/A2) * L2);
L_SE2 = @(t) L1;


% Full Flux (DO NOT CHANGE)
L_full = @(t) L_total;


% Graph Generation (DO NOT CHANGE)
hold on;
figure(1);
grid on;

title({"m_{1} = " + m1 + " M_{☉}, r_{1} = " + r1 + " R_{☉}, L_{1} = " + L1 + " L_{☉}", "m_{2} = " + m2 + " M_{☉}, r_{2} = " + r2 + " R_{☉}, L_{2} = " + L2 + " L_{☉}", "P = " + P/(24*60^2) + " d, a = " + sma/AU + " AU, i = " + i + "°"});
xlabel({"Seconds"})
ylabel("Solar Luminosities");
xlim([0, P]);
ylim([0, 1.2*L_total]);

% Piecewise Version (DO NOT CHANGE)
syms L_PE(t);
L_PE(t) = piecewise((t1 <= t & t <= t_total), L_PE1);
fplot(L_PE, "red");


% Piecewise Version (DO NOT CHANGE)
syms L_FF(t);
L_FF(t) = piecewise((t_total <= t & t <= P/2), L_full, (P/2 + t_total <= t & t <= P), L_full);
fplot(L_FF, "black");


% Piecewise Version (DO NOT CHANGE)
syms L_SE(t)
L_SE(t) = piecewise((P/2 <= t & t <= P/2 + t_total), L_SE1);
fplot(L_SE, "blue")


if b_h <= (r1 - r2)
    val = sqrt((((r1 - r2)^2 - b_h^2) * w^2 * R_Sol^2)/v_rel^2);
    t2 = (-asin(val) + phi)/w
    t3 = (asin(val) + phi)/w
    
    % Primary Eclipse Bottom Flux
    
    fplot(L_PE2, [t2, t3], "red");
    fplot(L_SE2, [P/2 + t2, P/2 + t3], "blue");
   
end

legend('Primary Eclipse', 'Full Flux', 'Secondary Eclipse', '', '');


hold off;
saveas(gcf, "binarycurve.png");
% --------------------------------------------
