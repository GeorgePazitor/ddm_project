clear all 
close all
clc
addpath('utils')
%% initialisation of parameters
L = 100000;   % (30) in mm
S = 1;     % (1)  mm^2
Fd = 10;   % (10) force in newton 
E = 2e5;   % (2e5)young's module in MPa

H = 1000;    % (10) H = lenght of a substructure (must be a divisor of L)
h = 10;     % (5)  h = lenght of an element inside a single substructure (must be a divisor of H)

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

%% kernel of Sp_s = rigid body modes of each substructure

Rb_l = cell(1,N);
for s = 1:N
    Sp_l{s};
    % null or eig always return a normalized vector 
    % null computes the kernel of the Sp(s) matrix -> rigid.b modes of s
    Rb_l{s} = null(Sp_l{s}, 'rational'); % not normalized rigid body modes 
    s;
    Rb_s = Rb_l{s};
    Rb_s;
end

%% computing Sd

Sd_l = cell(1, N);
for s=1:N
    Sd_l{s} = pinv(Sp_l{s}); 
end

Sd_diam = blkdiag(Sd_l{:});

bp_diam = cat(1, bp_l{:});

Sd = Ab_diam * Sd_diam * Ab_diam'; 

%% dual shur complement operators

Rb_diam = blkdiag(Rb_l{:});

G = Ab_diam * Rb_diam;  

e_diam = Rb_diam' * bp_diam;

bd = Ab_diam * Sd_diam * bp_diam; 


%% assemble all the necessary operators 
bp_diam = cat(1, bp_l{:});

Sp_diam = blkdiag(Sp_l{:});

Sd_diam = blkdiag(Sd_l{:}); %Sp plus diam

Rb_diam = blkdiag(Rb_l{:});

bp = A_diam * bp_diam;

I = cell(1, N);
for s = 1:N
    if s==1
        I{s} = eye(1);
    elseif s==N
        I{s} = eye(1);
    else
        I{s} = eye(2); % identity matrix for internal degrees of freedom
    end
end

M_diam = blkdiag(I{:}); % 1d bar: n elements -> n+1 dofs same as M 

Ab_tild  = (Ab_diam * inv(M_diam) * Ab_diam')' *  Ab_diam * inv(M_diam) ;

%Sd_tild = Ab_tild * Sd_diam * Ab_tild';

%Sd = Ab_diam * Sd_diam * Ab_diam';

%% Preconditioned conjugate gradient with Neumann preconditioner 
%initialization 

m=100;
epsilon = 1e-10;

beta_i = cell(1,m);
p = cell(1, m);
d = cell(1, m);

p_l = cell(1, N);
d_l = cell(1, N);
zi_l = cell(1, N);
ri_l = cell(1,N);
tmp_l = cell(1,N);
lambda_l = cell(1, N);

Q = eye(size(G,1));
P = eye(size(Q)) - Q * G * inv(G' * Q * G) * G';
lambda_i = - Q * G * inv(G' * Q * G) * e_diam;

% Distribute initialization of lambda to
% solve the local Neumann problem
lambda_diam = Ab_diam'*lambda_i;
lambda_l = vert_diam_to_s(lambda_diam, N);

for s=1:N
    ri_l{s} =  Sd_l{s} * lambda_l{s}; 
end

ri_diam = cat(1, ri_l{:});
ri = P'*(-bd-(Ab_diam*ri_diam));

% Distribute initialization of lambda to
% solve for the local Dirichlet preconditioner
ri_diam = Ab_tild'*ri;
ri_l = vert_diam_to_s(ri_diam, N);

for s=1:N
    tmp_l{s} = Sp_l{s}*ri_l{s};
end

tmp_diam = cat(1, tmp_l{:});
zi = P'* (Ab_tild*tmp_diam); %apply the projector to search on the good search space 

d{1}=zi;   %initialize the search direction

r0_norm = norm(ri);
if norm(ri)/r0_norm > epsilon 
    for i=1:m
        i
        % Distribute initialization of d to
        % solve the local Neumann problem
        d_diam = Ab_diam' * d{i};
        d_l = vert_diam_to_s(d_diam, N );
    
        for s=1:N     
            p_l{s} = Sd_l{s} * d_l{s};
        end
    
        p_diam= cat(1, p_l{:});
        p{i}= P' * (Ab_diam * p_diam); %apply the projector
    
        % global conjugate gradient on assembled residue
        alpha_i = (ri'*d{i}) / (d{i}'*p{i});  
    
        lambda_i = lambda_i + alpha_i * d{i};
        
        ri = ri - alpha_i * p{i};
    
        % Distribute the new residue to
        % solve for the local Dirichlet preconditioner
        ri_diam = Ab_tild' * ri;
        r_l = vert_diam_to_s(ri_diam, N );
    
        for s=1:N
            zi_l{s} = Sp_l{s}*r_l{s};
        end
    
        zi_diam = cat(1, zi_l{:});
        zi = P * (Ab_tild * zi_diam);
        
        % reorthogonalization of the search directions
        % necessary due to numerical incaccuracies
        update = zeros(size(zi));
        for j=1:i
            beta_i{j} = -(zi)'*p{j} / (d{j}'*p{j});
            update = update + beta_i{j}*d{j};     
        end
    
        d{i+1} = zi + update; %update the search direction
    
        %%%------------
        if norm(ri)/r0_norm < epsilon 
            break;
        end
        %%%------------
    end
    
end

%post processing - TODO, distribute the process 
alpha_diam = (G' * Q * G ) \ G' * Q * (- bd -  Sd * lambda_i );
ub_diam = Sd_diam * (bp_diam + Ab_diam' * lambda_i) + Rb_diam * alpha_diam;

ub = 1/2 * A_diam * ub_diam; % 1/multiplicity comes from A_diam*A_diam' = 1 / multiplicity* I (identity)
ub


