clear all 
close all
clc
addpath('utils');
%% initialisation of parameters
L = 100000;   % (30) in mm
S = 1;     % (1)  mm^2
Fd = 10;   % (10) force in newton 
E = 2e5;   % (2e5)young's module in MPa

H = 10000;    % (10) H = lenght of a substructure (must be a divisor of L)
h = 1000;     % (5)  h = lenght of an element inside a single substructure (must be a divisor of H)

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

%% Generate list containing the global stiffness matrices for each s, FEM
[K_l, f_l] = global_k_f_lists(N, H, E, S, Fd, n, node); 

%% Compute the internal,interface and combined submatrices and subvectors
[Kii_l, Kib_l, Kbi_l, Kbb_l, fi_l, fb_l] = internal_interface_partition(K_l, f_l, i_list, b_list, N);

%% Compute the local dirichlet operator Sp (and bp) for each substructure
% (Primal shur complement)

Sp_l = cell(1, N);
bp_l = cell(1, N);

for s=1:N    
    Sp_l{s} = Kbb_l{s}-Kbi_l{s}* (Kii_l{s} \ Kib_l{s});
    bp_l{s} = fb_l{s}-Kbi_l{s}* (Kii_l{s} \fi_l{s});
end

Rb_l = cell(1,N);
for s = 1:N
    Rb_l{s} = null(Sp_l{s}, 1e-6); % normalized version of the rigid body modes     
end

% use {:} to expand the list in its individual items
Sp_diam = blkdiag(Sp_l{:}); % diagonal concatenation
bp_diam = cat(1,bp_l{:});   % 1 for concatenation in first dim : vertical

%% ############## TODO ############# 
% build the discretized trace operator for each substructure and store them 
% into a list/array that you can initialize with the command cell
% as you can see in the following A and A bar cases

Sp = A_diam*Sp_diam*A_diam';

bp = A_diam*bp_diam;

% solution of primal shur complement with direct method
tic
ub = Sp\bp
toc

ub_diam=A_diam'*ub;

%% compute ui from ub 

ub_l=vert_diam_to_s(ub_diam, N);

ui_l = cell(1,N);
for s=1:N
    if s == 1        
        Kii = K_l{s}(i_list{s}(2:end), i_list{s}(2:end));
        Kib = K_l{s}(i_list{s}(2:end), b_list{s});
        fi = f_l{s}(i_list{s}(2:end));

        ui_l{s} = Kii\(fi-Kib*ub_l{s});

    else
        Kii = K_l{s}(i_list{s}, i_list{s});
        Kib = K_l{s}(i_list{s}, b_list{s});
        fi = f_l{s}(i_list{s});

        ui_l{s} = Kii\(fi-Kib*ub_l{s});
    end
    %ui_l{s};
end


%% kernel of Sp_s = rigid body modes of each substructure

R_l = cell(1,N);
for s = 1:N
    % null or eig always return a normalized vector 
    % null computes the kernel of the Sp(s) matrix -> rigid.b modes of s
    R_l{s} = null(Sp_l{s}); % normalized version of the rigid body modes 
    s;
    R_s = R_l{s};
    R_s;
end






