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

%% build discretized trace operator and assebly operators A and A bar

% TODO build discretized trace operator

A_diam = A_op(N);

Ab_diam = A_bar_op(N);

%% FEM for each substructure

node = [1:n        % matrix of n column vectors, each of which has the 
        2:n+1];    % corresponds to an element an its constitutive nodes  
%% Generate internal and interface node lists for each substructure
[i_list, b_list] = node_lists(N,n); 

%% Generate list containing the global stiffness matrices for each s
K_l = cell(1,N); 
f_l = cell(1,N);

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
%% compute the internal,interface and combined submatrices and subvectors

inv_Kii_l = cell(1, N);
Kii_l = cell(1, N);
Kib_l = cell(1, N);
Kbi_l = cell(1, N);
Kbb_l = cell(1, N);

fi_l = cell(1, N);
fb_l = cell(1, N);

for s=1:N    
    if s == 1 %elimination of the degree of freedom to impose ud
        %just exclude the first line and/or column for internal dofs 
        Kii_l{s} = K_l{s}(i_list{s}(2:end), i_list{s}(2:end));
        Kib_l{s}= K_l{s}(i_list{s}(2:end), b_list{s});
        Kbi_l{s} = K_l{s}(b_list{s}, i_list{s}(2:end));
        Kbb_l{s} = K_l{s}(b_list{s}, b_list{s});

        fi_l{s} = f_l{s}(i_list{s}(2:end));
        fb_l{s} = f_l{s}(b_list{s});

        inv_Kii_l{s} = inv(Kii_l{s});
    else

        Kii_l{s} = K_l{s}(i_list{s}, i_list{s});
        Kib_l{s} = K_l{s}(i_list{s}, b_list{s});
        Kbi_l{s} = K_l{s}(b_list{s}, i_list{s});
        Kbb_l{s} = K_l{s}(b_list{s}, b_list{s});
    
        fi_l{s} = f_l{s}(i_list{s});
        fb_l{s} = f_l{s}(b_list{s});

        inv_Kii_l{s} = inv(Kii_l{s});
    end
end

%% Distributed conjugate gradient to solve primal interface problem

%Initialization phase

ub = zeros(N-1, 1); %% n rows = n of boundary nodes
ub_diam = Ab_diam'*ub;

ub_l=vert_diam_to_s(ub_diam, N);

%solve the dirichlet problem (DOUBT ON HOW TO IMPOSE ud)
ui_l=cell(1,N);

rb_l=cell(1,N);
for s=1:N
    
        ui_l{s} = inv_Kii_l{s} * (fi_l{s} - Kib_l{s}*ub_l{s});
        
        Kbi_bb = cat(2, Kbi_l{s}, Kbb_l{s});   %OSS:possible optimisation is to avoid cat
        
        uib = cat(1, ui_l{s}, ub_l{s});
        
        rb_l{s} = -( Kbi_bb * uib - fb_l{s});
    
end

%global residue computation

rb_diam = cat(1,rb_l{:});
rb = A_diam*rb_diam;
db = rb;

db_l=cell(1,N);
di_l=cell(1,N);
Sp_db_l=cell(1,N);
m = 100;
for k=0:m
    k
    db_diam = A_diam'*db;
    db_l = vert_diam_to_s(db_diam, N );
    % solve local Dirichlet problem
    for s=1:N
        
        di_l{s} = -inv_Kii_l{s}*Kib_l{s}*db_l{s};
        
        %compute local matrix vector product
        Kbi_bb = cat(2, Kbi_l{s}, Kbb_l{s}); %OSS: possible optimisation is to avoid cat
        dib = cat(1, di_l{s}, db_l{s});
        
        Sp_db_l{s} = Kbi_bb*dib;
        s;
        Sp_db_l{s};
        
    end
    
    % Conjugate gradient update step

    Sp_db_diam = cat(1, Sp_db_l{:});
    Sp_db = A_diam*Sp_db_diam;             %Compute global matrix-vector product
    
    alpha = (rb'*rb) / (db'*Sp_db);        %Compute optimal step
   
    ub_next = ub + alpha*db;               %Compute iterate

    rb_next = rb - alpha*Sp_db;            %Compute residual
    
    beta =  (rb_next'*rb_next) / (rb'*rb); %Compute orthogonalization parameter 
    
    db = rb_next + beta*db;                %Update search direction

    rb = rb_next;
    ub = ub_next;
    %%%------------
    if norm(db) < 1e-6 
        break;
    end
    %%%------------
    
end


