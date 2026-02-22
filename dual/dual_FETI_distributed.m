clear all 
clc

%% initialisation of parameters

L = 1000;   % (30) in mm %k0 = 40 000
S = 1;     % (1)  mm^2
Fd = 10;   % (10) force in newton 
E = 2e5;   % (2e5)young's module in MPa

H = 100;    % (10) H = lenght of a substructure (must be a divisor of L)
h = 10;     % (5)  h = lenght of an element inside a single substructure (must be a divisor of H)

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
    %%% ---FOR THE SAKE OF ITERATIVE METHODS VERIFICATION---
    %if s == 2 
    %    f_l{s}(2) = Fd; % apply external force 
    %end
    %%% ----------------------------------------------------
end
%% compute the internal,interface and combined submatrices and subvectors

inv_Kii_l = cell(1, N);
Kii_l = cell(1, N);
Kib_l = cell(1, N);
Kbi_l = cell(1, N);
Kbb_l = cell(1, N);

fi_l = cell(1, N);
fb_l = cell(1, N);

Sp_l = cell(1, N);
bp_l = cell(1, N);


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
    %%% ---FOR THE SAKE OF ITERATIVE METHODS VERIFICATION---
    %elseif s == N %elimination of the degree to impose ud = 0 at the end
    %    Kii_l{s} = K_l{s}(i_list{s}(1:end-1), i_list{s}(1:end-1));
    %    Kib_l{s}= K_l{s}(i_list{s}(1:end-1), b_list{s});
    %    Kbi_l{s} = K_l{s}(b_list{s}, i_list{s}(1:end-1));
    %    Kbb_l{s} = K_l{s}(b_list{s}, b_list{s});

    %    fi_l{s} = f_l{s}(i_list{s}(1:end-1));
    %    fb_l{s} = f_l{s}(b_list{s});

    %    inv_Kii_l{s} = inv(Kii_l{s});
    %%% ----------------------------------------------------
    else

        Kii_l{s} = K_l{s}(i_list{s}, i_list{s});
        Kib_l{s} = K_l{s}(i_list{s}, b_list{s});
        Kbi_l{s} = K_l{s}(b_list{s}, i_list{s});
        Kbb_l{s} = K_l{s}(b_list{s}, b_list{s});
    
        fi_l{s} = f_l{s}(i_list{s});
        fb_l{s} = f_l{s}(b_list{s});

        inv_Kii_l{s} = inv(Kii_l{s});
    end

    Sp_l{s} = Kbb_l{s}- Kbi_l{s}*inv_Kii_l{s}*Kib_l{s};

    bp_l{s} = fb_l{s} - Kbi_l{s} * inv_Kii_l{s}* fi_l{s};
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


