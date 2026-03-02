clear all 
close all
clc
addpath('utils')
%% initialisation of parameters
L = 1000000;   % (30) in mm
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
Sd_l = cell(1, N);
bp_l = cell(1, N);

for s=1:N    

    Sp_l{s} = Kbb_l{s}-Kbi_l{s}*(Kii_l{s}\Kib_l{s});

    bp_l{s} = fb_l{s}-Kbi_l{s}*(Kii_l{s}\fi_l{s});

    Sd_l{s} = pinv(Sp_l{s});
end

%% kernel of Sp_s = rigid body modes of each substructure

Rb_l = cell(1,N);
for s = 1:N
    Rb_l{s} = null(Sp_l{s}); % normalized version of the rigid body modes     
end

%% assemble all the necessary operatiors 
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

M_diam  = blkdiag(I{:}); % 1d bar: n elements -> n+1 dofs 

A_tild  = (A_diam * M_diam * A_diam') \ A_diam * M_diam;

Sp_tild = A_tild * Sd_diam * A_tild';

Sp = A_diam * Sp_diam * A_diam';

G_tilde = A_tild * Rb_diam; 

%% conjugate gradient 
%initialization 

m=3*N;
epsilon = 1e-10;
beta_i = cell(1,m);
p = cell(1,m);
d = cell(1,m);

% here if I compute the inverse with \ I wil get the wrong solution
ui = zeros(size(bp)); %pinv(Sp)* bp;

ri = (bp - Sp * ui ); % = P'*bp

d{1} = ri;

residuals = zeros(1, m);

r0_norm = norm(ri);
if r0_norm > epsilon

    for i=1:m 
        i;
        p{i} = Sp*d{i};
    
        alpha_i = (ri'*d{i}) / (d{i}'*p{i});  %compute the optimal step
    
        ui_next = ui + alpha_i * d{i};
    
        ri_next = ri - alpha_i * p{i};
        
        %update = zeros(size(zi));
    
        %for j=1:i
        %    beta_i{j} = -(ri_next)'*p{j} / (d{j}'*p{j});
        %    update = update + beta_i{j}*d{j}; 
        %end
    
        beta_i = -ri_next'*p{i} / (d{i}'*p{i});

        d{i+1} = ri_next + beta_i*d{i};%update; %update the search direction
    
        ri = ri_next;
        ui = ui_next;

        residuals(i) = norm(ri)/r0_norm;
        %%%--------convergence criterion-----
        if norm(ri)/r0_norm < epsilon 
            break;
        end
        %%%----------------------------------
    end
end
ub = ui;
i;
ub;


%% Precond conjugate gradient with neumann preconditioner
%initialization 

m=3*N;
epsilon = 1e-10;
beta_i = cell(1,m);
p = cell(1,m);
d = cell(1,m);

% here if I compute the inverse with \ I wil get the wrong solution
ui = zeros(size(bp)); %pinv(Sp)* bp;

ri = (bp - Sp * ui ); % = P'*bp
d{1} = ri;
zi = Sp_tild *ri;
d{1} = zi;

residuals_prec = zeros(1, m);
r0_norm = norm(ri);
if r0_norm > epsilon

    for i=1:m 
        i;
        p{i} = Sp*d{i};
    
        alpha_i = (ri'*d{i}) / (d{i}'*p{i});  %compute the optimal step
    
        ui_next = ui + alpha_i * d{i};
    
        ri_next = ri - alpha_i * p{i};
    
        zi_next = Sp_tild * ri_next;
    
            
        beta_i = -zi_next'*p{i} / (d{i}'*p{i});

        d{i+1} = zi_next + beta_i * d{i}; %update the search direction
    
        ri = ri_next;
        ui = ui_next;

        residuals_prec(i) = norm(ri)/r0_norm;
        %%%--------convergence criterion-----
        if norm(ri)/r0_norm < epsilon 
            break;
        end
        %%%----------------------------------
    end
end
ub = ui;
i;
ub;

plot(1:m, residuals, Color="red");
hold on
plot(1:m, residuals_prec, Color="blue");
legend("CG", "CG + perconditioner")
xlabel("Iteration i");
ylabel("norm\_ri / norm\_r0");
title("Relative residual comparison");


