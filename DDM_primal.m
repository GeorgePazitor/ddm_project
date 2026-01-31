clear all 
clc

%% initialisation of parameters

L = 30;   % in mm
S = 1;    % mm^2
Fd = 10;   % force in newton 
E = 2e5;   % young's module in MPa

H = 10;    % H = lenght of a substructure (must be a divisor of L)
h = 5;     % h = lenght of an element inside a single substructure (must be a divisor of H)


n = H/h; % number of element per substructure
N = L/H; % number of substructures 

%% Analytical solution
x = linspace(0, L, 1000); %define a discretized domain for x
u_x = @(x) (Fd/(E*S))*x ; %define the analytical solution in displacement 
stress_xx = Fd/S * ones(1000, 1); %define the analytical solution in stress

%figure(1) 
%plot(x, u_x(x), 'g-' )
%hold on

%figure(2)
%plot(x, stress_xx, 'b-')
%hold on

%% FEM for each substructure
%F = zeros(n+1, 1); % initialisation of the shape of the force vector
%F(n+1) = Fd;       % external force applied to the last node 



node = [1:n        % matrix of n column vectors, each of which has the 
        2:n+1];    % corresponds to an element an its constitutive nodes  

i_list = cell(1,N); %internal nodes 
b_list = cell(1,N); %interface nodes

%generation of internal and interface node lists
for s=1:N % for each substructure
    for e=1:n+1          % for each node in substructure
        if e == 1
            if s == 1 % first node of first substructure is internal 
                i_list{s} = [i_list{s} e];
            else
                b_list{s} = [b_list{s} e]; %otherwise it's interface
            end
        elseif e == n+1
            if s == N % last node of first substructure is internal 
                i_list{s} = [i_list{s} e];
            else
                b_list{s} = [b_list{s} e]; %otherwise it's interface
            end
        else
            i_list{s} = [i_list{s} e]; %otherwise it's an internal node
        end
    end 
end

K_l = cell(1,N); 
f_l = cell(1,N);
%generation of the list containing the global stiffness matrices for each s
for s=1:N % for each substructure s

    x = (H*(s-1)):H/n:(H*s);   % redefine x domain: 
    K = zeros(n+1);    % initialize the global stiffness matrix

    for e=1:n          % for each element 

        le = x(node(2, e)) - x(node(1,e));    % element lenght x2-x1
        ke = E*S/le * (eye(2) - flip(eye(2)));% local stiffness matrix

        % assemble global stiffness matrix with local contribution of ke
        K(node(1, e), node(1,e)) =  K(node(1, e), node(1,e)) + ke(1,1);
        K(node(1, e), node(2,e)) =  K(node(1, e), node(2,e)) + ke(1,2);
        K(node(2, e), node(1,e)) =  K(node(2, e), node(1,e)) + ke(2,1);
        K(node(2, e), node(2,e)) =  K(node(2, e), node(2,e)) + ke(2,2);
    end
    K_l{s} = K; % store the global stiffness matrix for each substructure
    f_l{s} = zeros(n+1, 1);
    if s == N 
        f_l{s}(end) = Fd; % apply external force to the last node of the last substructure
    end
end

Sp_l = cell(1, N);
bp_l = cell(1, N);
% find the local dirichlet operator Sp (and bp) for each substructure
for s=1:N
    if s == 1 %elimination of the degree of freedom to impose ud
        %just exclude the first line and column for internal dofs 
        Kii = K_l{s}(i_list{s}(2:end), i_list{s}(2:end));
        Kib = K_l{s}(i_list{s}(2:end), b_list{s});
        Kbi = K_l{s}(b_list{s}, i_list{s}(2:end));
        Kbb = K_l{s}(b_list{s}, b_list{s});

        fi = f_l{s}(i_list{s}(2:end));
        fb = f_l{s}(b_list{s});
    else

        Kii = K_l{s}(i_list{s}, i_list{s});
        Kib = K_l{s}(i_list{s}, b_list{s});
        Kbi = K_l{s}(b_list{s}, i_list{s});
        Kbb = K_l{s}(b_list{s}, b_list{s});
    
        fi = f_l{s}(i_list{s});
        fb = f_l{s}(b_list{s});
    end
    inv_Kii = inv(Kii);
    Sp_l{s} = Kbb-Kbi*inv_Kii*Kib;

    bp_l{s} = fb-Kbi*inv_Kii*fi;

end

% use {:} to expand the list in its individual items
Sp_diam = blkdiag(Sp_l{:}); % diagonal concatenation
bp_diam = cat(1,bp_l{:}); % 1 for concatenation in first dim : vertical

%% ############## TODO ############# 
% build the discretized trace operator for each substructure and store them 
% into a list/array that you can initialize with the command cell
% as you can see in the following A and A bar cases

%%
% build A diamond
A_l = cell(1,N);        
for s=1:N
    %for the first and the last substruct it's just a vector A(s) 
    % has 1 column (1 single interaction per substructure)
    if s == 1 || s == N   
        A_l{s} = zeros(N-1, 1); 
    %for the first and the last substruct it's just a vector A(s) has 2 col
    else                 
        A_l{s} = zeros(N-1, 2);
    end
end

for c=1:(N-1)
    if c == 1
        A_l{c}(c, 1) = 1;
        A_l{c+1}(c, 1) = 1; 
    else
        A_l{c}(c, 2) = 1; 
        A_l{c+1}(c, 1) = 1; 
    end
    
end

A_diam = cat(2, A_l{:}) 

% build A bar diamond
Ab_l = cell(1,N);
for s=1:N
    if s == 1 || s == N
        Ab_l{s} = zeros(N-1, 1); 
    else
        Ab_l{s} = zeros(N-1, 2);
    end
end

for c=1:(N-1)
    if c == 1
        Ab_l{c}(c, 1) = 1;
        Ab_l{c+1}(c, 1) = -1; 
    else
        Ab_l{c}(c, 2) = 1; 
        Ab_l{c+1}(c, 1) = -1; 
    end
    
end

Ab_diam = cat(2, Ab_l{:}) % 2: horizontal concatenation

Sp = A_diam*Sp_diam*A_diam'
bp = A_diam*bp_diam

% solution of primal shur complement with direct method
ub = Sp\bp

ub_diam=A_diam'*ub;

ui_l = cell(1:N);
for s=1:N
    if s == 1
        ub_s = ub_diam(s, 1);
        Kii = K_l{s}(i_list{s}(2:end), i_list{s}(2:end));
        Kib = K_l{s}(i_list{s}(2:end), b_list{s});
        fi = f_l{s}(i_list{s}(2:end));

        ui_l{s} = inv(Kii)*(fi-Kib*ub_s);

    else
        if s==N
            ub_s = ub_diam((s*2)-2, 1);
        else
            ub_s = ub_diam(((s*2)-2):((s*2)-1), 1);
        end
        Kii = K_l{s}(i_list{s}, i_list{s});
        Kib = K_l{s}(i_list{s}, b_list{s});
        fi = f_l{s}(i_list{s});

        ui_l{s} = inv(Kii)*(fi-Kib*ub_s);
    end
    ui_l{s};
end

ui_diam = cat(1, ui_l{:});
ui = A_diam*ui_diam;
%u = [ui; ub];       

% kernel of Sp_s = rigid body modes of each substructure

R_l = cell(1,N);
for s = 1:N
    % null or eig always return a normalized vector 
    % null computes the kernel of the Sp(s) matrix -> rigid.b modes of s
    R_l{s} = null(Sp_l{s}); % normalized version of the rigid body modes 
    %R_l{s}
end






