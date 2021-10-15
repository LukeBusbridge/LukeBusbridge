%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%ROAD SIGN VERIFICATION%%%%
%%%%%%LUKE BUSBRIDGE 2021%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Pole Parameters
radius = 200; %radius of pole in mm
thickness = 10; %thickness of the pole (mm
area = 2*pi*radius*thickness; %calculated cross sectional area of the pole
I = pi*radius^3*thickness; %second momnet of area (mm^4)
L_v = 13300/3 %vertical length of the pole (mm)
L_h = 3900 %horizontal length of the pole (mm)
E = 205e3; %modulus of elasticity of the pole (MPa)
density_pole = 7870e-12; %density of pole material (t/mm^4)
density_sign = 2680e-12; %density of sign material (t/mm^4)
mass_man = 100; % mass of man (kg)

%DISPLACEMENT RESPONSE

%Calculation of forces and moments
total_length = L_v+L_h %length of the entire pole
%Self weight of pole and calculated distributed load
w_pole = L_h*area*density_pole*9810; %N
dist_pole = w_pole/L_h; %(N/mm)
F_pole = -dist_pole*3900/2; %determine the applied force from pole weight
gravity_moment_pole = -dist_pole*3900^2/12; %moment from pole weight
%self weight of the sign and calculated distributed load
w_sign = 4.93e6*1.6*density_sign*9810; %N
dist_sign = w_sign/2900; %(N/mm)
F_sign = -dist_sign*2900/2; %determine the applied force from sign weight
gravity_moment_sign = -dist_sign*2900^2/12; %moment from sign weight
W_man = mass_man*1e-3*-9810; %weight of man on the end of pole

BeamVert = zeros(9,9); %create an array of zeros for the verical beam
%Assign the beam stiffness matrix
BeamVert([1:6],[1:6]) = [E*area/L_v, 0, 0, -E*area/L_v, 0, 0;...
        0, 12*E*I/L_v^3, 6*E*I/L_v^2, 0, -12*E*I/L_v^3, 6*E*I/L_v^2;...
        0, 6*E*I/L_v^2, 4*E*I/L_v, 0, -6*E*I/L_v^2, 2*E*I/L_v;...
        -E*area/L_v, 0, 0, E*area/L_v, 0, 0;...
        0, -12*E*I/L_v^3, -6*E*I/L_v^2, 0, 12*E*I/L_v^3, -6*E*I/L_v^2;...
        0, 6*E*I/L_v^2, 2*E*I/L_v, 0, -6*E*I/L_v^2, 4*E*I/L_v;];
 
BeamHorizontal = zeros(9,9); %create an array of zeros for the horizontal beam
%Assign the beam stiffness matrix
BeamHorizontal([4:9],[4:9]) = [E*area/L_h, 0, 0, -E*area/L_h, 0, 0;...
        0, 12*E*I/L_h^3, 6*E*I/L_h^2, 0, -12*E*I/L_h^3, 6*E*I/L_h^2;...
        0, 6*E*I/L_h^2, 4*E*I/L_h, 0, -6*E*I/L_h^2, 2*E*I/L_h;...
        -E*area/L_h, 0, 0, E*area/L_h, 0, 0;...
        0, -12*E*I/L_h^3, -6*E*I/L_h^2, 0, 12*E*I/L_h^3, -6*E*I/L_h^2;...
        0, 6*E*I/L_h^2, 2*E*I/L_h, 0, -6*E*I/L_h^2, 4*E*I/L_h];    

%Transpose the verical beam Matrix
alpha= pi/2; %angle of first beam (rad)
%Calculate transpose Matrix
T= zeros(9,9);
T([1:6],[1:6])=[cos(alpha), sin(alpha), 0, 0, 0, 0;...
    -sin(alpha), cos(alpha), 0, 0, 0, 0;...
    0, 0, 1, 0, 0, 0;...
    0, 0, 0, cos(alpha), sin(alpha), 0;
    0, 0, 0, -sin(alpha), cos(alpha), 0;
    0, 0, 0, 0, 0, 1];
   
Vert= T*BeamVert; %perform transpose to veritcal
Total = Vert+BeamHorizontal; %combine the verical and horizontal matrices
T_reduced = Total([4:9],[4:9]); %reduce the stiffness matrix from boundary conditions

%Allocating correct forces in the force matrix [fx2, fy2, M2, fx3, fy3, M3] 
Forces = [F_pole+F_sign;0;-gravity_moment_pole-gravity_moment_sign;0;...
    F_pole+F_sign+W_man; gravity_moment_pole+gravity_moment_sign];
response = linsolve(T_reduced,Forces); %solve for the dusplacements [u2,v2,thi2,u3,v3,thi3]

VerticalDeflection = response(5) %show the vertical deflection at the end of the pole

%VIBRATIONAL RESPONSE
syms l %defile the l variable
eqn = (12-156*l)*(4*total_length^2-4*total_length^2*l)-(-6*total_length+22*l*total_length)*(-6*total_length+22*total_length*l) == 0; %set out the equation;
lambda = solve(eqn,l); %solve for lambda;
fq1 = sqrt(420*lambda(1))/(2*pi)*(sqrt(E*I/(density_pole*area*total_length^4))) %find the first natural frequency
