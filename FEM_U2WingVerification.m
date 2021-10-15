%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%BENDING U1 WING VERIFICATION%%%%
%%%%%%%%%LUKE BUSBRIDGE 2021%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sympref('FloatingPointOutput',true)
clear all;
%Apply Lift to wing?
prompt = 'Do you want to apply Lift? Y/N: ';
apply_lift = input(prompt,'s');
%material selection
prompt = 'aluminium or titanium? A/T: ';
material = input(prompt,'s');

%Wing Parameters
chord=3; %meters
airfoil_max_thickness = 0.065; %0.121 for rae2822, 0.065 for raf15
height=airfoil_max_thickness*chord; %meters
thickness=0.050; %meters  0.003 for aluminium, 0.0025 for titanium
span = 14; %meters
fuel_loc = 4; % m from inboard
fuel_mass = 10887/2; %kg
aircraft_mass=7257/2; %weight supported by one wing
gravity = 9.81; %m/s^2
percent_rectangle = 0.2; %percent of the airfoil section to be considered a rectangle

%Check whether lift should be applied or not
if apply_lift== "Y"
    lift_factor=1;
    title_text = "with";
else
    lift_factor=0;
    title_text = "with no";
end


%Material Properties
if material == "A"
%http://www.matweb.com/search/DataSheet.aspx?MatGUID=da98aea5e9de44138a7d28782f60a836&ckck=1
    material_name = "aluminium";
    density = 2810; %kg/m^3
    poisson = 0.33;
    modulus_elasticity = 71.7*10^9; %pa
    yield = 103*10^6;
else 
%http://www.matweb.com/search/DataSheet.aspx?MatGUID=66a15d609a3f4c829cb6ad08f0dafc01
    material_name = "titanium";
    density = 4500; %kg/m^3
    poisson = 0.34;
    modulus_elasticity = 116*10^9; %pa
    yield = 140*10^6;
end
yield_FOS = yield/2;

%Wing Properties
percent_ellipse = 1-percent_rectangle;
a=chord+thickness; %set the chord as a
b=height+thickness; %set the height as be
a1=a-2*thickness; %internal chord
b1=b-2*thickness; %internal height
area_r=a*b-a1*b1; %calculate area of a hollow rectangle
area_e= pi/4*(a*b-a1*b1);
area = (area_r*percent_rectangle+area_e*percent_ellipse);
volume= area*span; %volume of the wing
wing_mass= volume*density; % find the mass of the wing
Ix_r= a*b^3/12-a1*b1^3/12; %second moment of area for thin rectangle
Ix_e= pi*(a*b^3-a1*b1^3)/64;
Ix = (Ix_r*percent_rectangle+Ix_e*percent_ellipse);

%Calculation of loads (up +, right +, AC +)
total_mass = aircraft_mass+fuel_mass+wing_mass; %total mass of 50% of the aircraft
lift_load = lift_factor*total_mass*gravity; %lift requred for each wing (Newtons)
distributed_lift = lift_load/span; %N/m
lift_force_wing = distributed_lift*span/2; %N
lift_moment_wing = distributed_lift*span^2/12; %direction needs checking %Nm
gravity_load = gravity*wing_mass; %N
gravity_load_distributed = gravity_load/span; %N/m
gravity_force_wing = -gravity_load_distributed*span/2; %N
gravity_moment_wing = -gravity_load_distributed*span^2/12; %Nm
fuel_load = -fuel_mass*gravity; %N

%Stiffness Matrx
E = modulus_elasticity; %assign youngs modulus
L1=fuel_loc; %inboard length
L2=span-fuel_loc; %outboard length
K1=E*Ix/L1^3; K2=E*Ix/L2^3; % define stiffness constants
%Calculation of stiffness matrix
k = [12*K1,6*L1*K1,-12*K1,6*L1*K1,0,0;...
  6*L1*K1,4*L1^2*K1,-6*L1*K1,2*L1^2*K1,0,0;...
  -12*K1,-6*L1*K1,12*K1+12*K2,-6*L1*K1+6*L2*K2,-12*K2,6*L2*K2;...
  6*L1*K1,2*L1^2*K1,-6*L1*K1+6*L2*K2,4*L1^2*K1+4*L2^2*K2,-6*L2*K2,2*L2^2*K2;...
  0,0,-12*K2,-6*L2*K2,12*K2,-6*L2*K2;...
  0,0,6*L2*K2,2*L2^2*K2,-6*L2*K2,4*L2^2*K2+4];

%Reduced stiffness matrix based on boundary conditions
reducedK = [k(15),k(16),k(17),k(18);...
    k(21),k(22),k(23),k(24);...
    k(27),k(28),k(29),k(30);...
    k(33),k(34),k(35),k(36)];

%Apply boundary conditions
v1 = 0; phi1 = 0;

%Known Forces and moments
syms fy1 M1
fy2 = fuel_load;
fy3 = lift_force_wing+gravity_force_wing;
M2 = 0;
M3 = gravity_moment_wing+lift_moment_wing;
f = [fy1;M1;fy2;M2;fy3;M3]; %creation of corce matrix
reducedf = [f(3);f(4);f(5);f(6)]; %reduced force matrix
U_reduced = linsolve(reducedK,reducedf); %solve foe displacements

%Assign the solved values
v2 = U_reduced(1); phi2= U_reduced(2); v3= U_reduced(3); phi3= U_reduced(4);

%Full U matrix
U = [v1; phi1 ; v2; phi2; v3; phi3];
Forces= k*U; %ccalculate reaction forces

%%%%Figure Creation%%%%

%Sketch the wing
clf('reset')
hold off
yyaxis left
x_coords = [0,fuel_loc,span]; y_coords = [U(1),U(3),U(5)]; %wing coords
plot(x_coords,y_coords,'DisplayName','Deflected Wing');
hold on;
xlim([-1 span+1]); ylim([-2 2]); fuse_rad = 0.75; %defining circle features
viscircles([-fuse_rad,0],fuse_rad,'Color','blue');
xL = xlim;
line(xL, [0 0],'Color','black','LineStyle','--','DisplayName','Horizontal'); 
line(xL, [-1 -1],'Color','red','LineStyle','-.','LineWidth',1,'DisplayName','Max Allowed Deflection');
legend; 
title("Wing " + title_text + " lift applied (" + material_name + ")");
xlabel('wing span (meters)') ;
ylabel('wing deflection (meters)') ;

%%%%Vibration Calculation%%%%
%considering the beam to have only two nodes, no mass
syms l %defile the l variable
eqn = (12-156*l)*(4*span^2-4*span^2*l)-(-6*span+22*l*span)*(-6*span+22*span*l) == 0; %set out the equation;
lambda = solve(eqn,l); %solve for lambda;
fq1 = sqrt(420*lambda(1))/(2*pi)*(sqrt(E*Ix/(density*area*span^4))) %find the first natural frequency
fq2 = sqrt(420*lambda(2))/(2*pi)*(sqrt(E*Ix/(density*area*span^4))) %find the second natural frequency
