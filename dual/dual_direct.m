clear all 
close all
clc
addpath('utils')
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
node = [1:n        % matrix of n column vectors, each of which has the 
        2:n+1];    % corresponds to an element an its constitutive nodes  
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

    Sp_l{s} = Kbb_l{s}-Kbi_l{s}*(Kii_l{s}\Kib_l{s});

    bp_l{s} = fb_l{s}-Kbi_l{s}*(Kii_l{s}\fi_l{s});
end

% use {:} to expand the list in its individual items
Sp_diam = blkdiag(Sp_l{:}); % diagonal concatenation
bp_diam = cat(1,bp_l{:});   % 1 for concatenation in first dim : vertical

Sp = A_diam*Sp_diam*A_diam';
bp = A_diam*bp_diam ;

%% kernel of Sp_s = rigid body modes of each substructure

Rb_l = cell(1,N);
for s = 1:N
    Sp_l{s};
    % null or eig always return a normalized vector 
    % null computes the kernel of the Sp(s) matrix -> rigid.b modes of s
    Rb_l{s} = null(Sp_l{s}, 1e-6); % not normalized rigid body modes 
    s;
    Rb_s = Rb_l{s};
    Rb_s;
end

%% computing Sd

Sd_l = cell(1, N);
for s=1:N
    Sd_l{s} = pinv(Sp_l{s},  1e-6); 
end

Sd_diam = blkdiag(Sd_l{:});

Sd = Ab_diam * Sd_diam * Ab_diam'; 

%% solution of dual shur complement with direct method

Rb_diam = blkdiag(Rb_l{:});

G = Ab_diam * Rb_diam; 

e_diam = Rb_diam' * bp_diam; 

bd = Ab_diam * Sd_diam * bp_diam; 

B = [Sd, G;
    G', zeros(size(G,2))];

a = cat(1, -bd, -e_diam);

y = B \ a;

lamb = y(1:size(G,1)); 
alphab_diam  =  y(size(G,1)+1:end); 

lamb_diam = Ab_diam' * lamb; 

ub_diam = Sd_diam * (bp_diam + lamb_diam) + Rb_diam * alphab_diam;

ub_diam;

ub = 1/2 * A_diam * ub_diam;
% 1/multiplicity comes from A_diam*A_diam' = 1 / multiplicity* I (identity)

ub



