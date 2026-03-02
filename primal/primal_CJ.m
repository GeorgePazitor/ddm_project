clear all 
close all
clc
addpath('utils');
%% initialisation of parameters
L = 1000;   % (30) in mm
S = 1;     % (1)  mm^2
Fd = 10;   % (10) force in newton 
E = 2e5;   % (2e5)young's module in MPa

H = 100;    % (10) H = lenght of a substructure (must be a divisor of L)
h = 10;     % (5)  h = lenght of an element inside a single substructure (must be a divisor of H)

n = H/h; % number of element per substructure
N = L/H; % number of substructures 
%% build assebly operators A and A bar

A_diam = A_op(N);

Ab_diam = A_bar_op(N);

%% FEM for each substructure
node = [1:n        % matrix of n column vectors, each of which 
        2:n+1];    % corresponds to an element and its constitutive nodes  
%% Generate internal and interface node lists for each substructure
[i_list, b_list] = node_lists(N,n); 

%% Generate list containing the global stiffness matrices for each s

[K_l, f_l] = global_k_f_lists(N, H, E, S, Fd, n, node); 

%% Compute the internal,interface and combined submatrices and subvectors

[Kii_l, Kib_l, Kbi_l, Kbb_l, fi_l, fb_l] = internal_interface_partition(K_l, f_l, i_list, b_list, N);

%% Distributed conjugate gradient to solve primal interface problem

%Initialization phase
m = N;
epsilon = 1e-10;
ub = zeros(N-1, 1); %% n rows = n of boundary nodes
ub_diam = Ab_diam'*ub;

ub_l=vert_diam_to_s(ub_diam, N);

%solve the dirichlet problem 
ui_l=cell(1,N);

residuals_dist = zeros(1, m);

rb_l=cell(1,N);
for s=1:N
    
        ui_l{s} = Kii_l{s} \ (fi_l{s} - Kib_l{s}*ub_l{s});
        
        Kbi_bb = cat(2, Kbi_l{s}, Kbb_l{s});   

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

r0_norm = norm(rb);
k = 0;
while 1
%for k=1:N+1
    k = k+1;
    k;
    db_diam = A_diam'*db;
    db_l = vert_diam_to_s(db_diam, N );
    % solve local Dirichlet problem
    for s=1:N
        
        di_l{s} = -Kii_l{s}\(Kib_l{s}*db_l{s});
        
        %compute local matrix vector product
        Kbi_bb = cat(2, Kbi_l{s}, Kbb_l{s});
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
    residuals_dist(k)=norm(rb)/r0_norm;
    %%%--------convergence criterion-----
    if norm(rb)/r0_norm < epsilon 
        break;
    end
    %%%----------------------------------
    
end
k;
ub


%plot(1:m, residuals_dist, Color="green");
%xlabel("Iteration i");
%ylabel("norm\_ri / norm\_r0");
%title("relative residual distributed CG");
