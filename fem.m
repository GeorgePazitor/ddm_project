clear all 
clc

%% Finite element method for 1D bar with linear displacement interpolation

L = 100;   % in mm
S = 10;    % mm^2
Fd = 10;   % force in newton 
E = 2e5;   % young's module in MPa

% for the finite element method
n = 4;    % n = number of elements -> n+1 nodes

x = linspace(0, L, 1000); %define a discretized domain for x

u_x = @(x) (Fd/(E*S))*x ; %define the analytical solution in displacement 

stress_xx = Fd/S * ones(1000, 1); %define the analytical solution in stress

figure(1) 
plot(x, u_x(x), 'g-' )
hold on

figure(2)
plot(x, stress_xx, 'b-')
hold on

F = zeros(n+1, 1); % initialisation of the shape of the force vector
F(n+1) = Fd;       % external force applied to the last node 

x = 0:100/n:100;   % redefine x domain: 

node = [1:n        % matrix of n column vectors, each of which 
        2:n+1];    % corresponds to an element and its nodes  
                   % es: first elem is [ 1 2 ]^T of node 1 and 2

K = zeros(n+1);    % initialize the global stiffness matrix

for e=1:n          % for each element 
    e
    le = x(node(2, e)) - x(node(1,e));    % element lenght x2-x1
    ke = E*S/le * (eye(2) - flip(eye(2)))% local stiffness matrix
    % assemble the global stiffness matrix with local contribution of ke
    K(node(1, e), node(1,e)) =  K(node(1, e), node(1,e)) + ke(1,1);
    K(node(1, e), node(2,e)) =  K(node(1, e), node(2,e)) + ke(1,2);
    K(node(2, e), node(1,e)) =  K(node(2, e), node(1,e)) + ke(2,1);
    K(node(2, e), node(2,e)) =  K(node(2, e), node(2,e)) + ke(2,2);

end
K

u_1 = 0; % boundary conditions

% to solve the linear system -> remoove dof where ud = 0 
F_1 = F(2:n+1);        
K_1 = K(2:n+1, 2:n+1); % elimination of first row/first column

u_int = K_1\F_1;       % to solve a linear system K*u=f we can multiply 
                       % by the inverse of K on both sides
                       % so we have now u=K^-1*f, the backslash is used in
                       % matlab to solve linear systems and is placed with 
                       % a similar logic as the inverse of a matrix  

u_tot = [u_1           % full displacement solution
         u_int];

stress = zeros(n, 1);
for e=1:n
   le = (x(node(2, e)) - x(node(1,e)));
   % stress formula for 1D axially loaded bar
   % the ' symbol after a vector means transposed
   stress(e) = E/le*[ -1 1 ]*[ u_tot(node(1, e)) u_tot(node(2, e))]';

end


figure(1)
plot(x, u_tot, 'ro')

legend('Analytical Displacement', 'FEM Displacement' );
xlabel('Position (mm)');
ylabel('Displacement (mm)');
title('Displacement Comparison');

figure(2)
plot(x(1:n), stress, 'ro')
legend('Analytical Stress', 'FEM Stress');
xlabel('Position (mm)');
ylabel('Stress (MPa)');
title('Stress Comparison');


