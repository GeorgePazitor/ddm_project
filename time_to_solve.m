clear all 
clc

L = 100;   % in mm
S = 10;    % mm^2
Fd = 10;   % force in newton 
E = 2e5;   % young's module in MPa

d = 1000;  %number of discretization points for the analytical solution

% for the finite element method
n_min_dofs = 1000;
n_max_dofs = 20000;    % n = number of elements -> n+1 nodes

time_taken_fem = zeros(1,n_max_dofs/n_min_dofs);
time_taken_primal = zeros(1,n_max_dofs/n_min_dofs);
time_taken_dist_primal_CG = zeros(1,n_max_dofs/n_min_dofs);
j=1;
for i=n_min_dofs:n_min_dofs:n_max_dofs
%% FEM -------------------------------------------------------------------
    F = zeros(i+1, 1); % initialisation of the shape of the force vector
    F(i+1) = Fd;       % external force applied to the last node 
    
    x = 0:L/i:L;   % redefine x domain: 
    
    node = [1:i        % matrix of n column vectors, each of which 
            2:i+1];    % corresponds to an element and its nodes  
                       % es: first elem is [ 1 2 ]^T of node 1 and 2
    
    K = zeros(i+1);    % initialize the global stiffness matrix
    
    for e=1:i          % for each element 
        e;
        le = x(node(2, e)) - x(node(1,e));    % element lenght x2-x1
        ke = E*S/le * (eye(2) - flip(eye(2)));% local stiffness matrix
        % assemble the global stiffness matrix with local contribution of ke
        K(node(1, e), node(1,e)) =  K(node(1, e), node(1,e)) + ke(1,1);
        K(node(1, e), node(2,e)) =  K(node(1, e), node(2,e)) + ke(1,2);
        K(node(2, e), node(1,e)) =  K(node(2, e), node(1,e)) + ke(2,1);
        K(node(2, e), node(2,e)) =  K(node(2, e), node(2,e)) + ke(2,2);
    
    end
    K;
    u_1 = 0; % boundary conditions
    % to solve the linear system -> remoove dof where ud = 0 
    F_1 = F(2:i+1);        
    K_1 = K(2:i+1, 2:i+1); % elimination of first row/first column
    t_start = tic;
    u_int = K_1\F_1;       % to solve a linear system K*u=f we can multiply 
    time_taken_fem(j) = toc(t_start) ;  % by the inverse of K on both sides

    % ---------------------------------------------------------------------
    %% Primal direct ------------------------------------------------------
   
    L = i;   % (30) in mm
    S = 1;     % (1)  mm^2
    Fd = 10;   % (10) force in newton 
    E = 2e5;   % (2e5)young's module in MPa
    
    H = i*0.2;    % (10) H = lenght of a substructure (must be a divisor of L)
    h = 1;     % (5)  h = lenght of an element inside a single substructure (must be a divisor of H)
    
    n = H/h; % number of element per substructure
    N = L/H; % number of substructures 
        
    node = [1:n        % matrix of n column vectors, each of which
            2:n+1];    % corresponds to an element and its constitutive nodes  
    
    % Generate internal and interface node lists for each substructure
    [i_list, b_list] = node_lists(N,n); 
    
    % Generate list containing the global stiffness matrices for each s, FEM
    K_l = cell(1,N); 
    f_l = cell(1,N);
    
    for s=1:N % for each substructure s
        x = (H*(s-1)):H/n:(H*s);   % redefine x domain:takes into account the 
                                   % position in the reference system 
    
        K = zeros(n+1);    % initialize the global stiffness matrix
    
        for e=1:n          % for each element 
    
            le = x(node(2, e)) - x(node(1,e));    % element lenght x2-x1
            ke = E*S/le * (eye(2) - flip(eye(2)));% local stiffness matrix
            % assemble global stiffness matrix with local contribution of ke
            K(node(1, e), node(1, e)) =  K(node(1, e), node(1,e)) + ke(1,1);
            K(node(1, e), node(2, e)) =  K(node(1, e), node(2,e)) + ke(1,2);
            K(node(2, e), node(1, e)) =  K(node(2, e), node(1,e)) + ke(2,1);
            K(node(2, e), node(2, e)) =  K(node(2, e), node(2,e)) + ke(2,2);
        end
        K_l{s} = K; % store the global stiffness matrix for each substructure
        f_l{s} = zeros(n+1, 1); %initialize the force vector for each s
        if s == N 
            f_l{s}(end) = Fd; % apply external force to the last node of the last substructure
        end
    end
    Sp_l = cell(1, N);
    bp_l = cell(1, N);
    
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
            t_start = tic;
            inv_Kii_l{s} = inv(Kii_l{s});
            Sp_l{s} = Kbb_l{s}-Kbi_l{s}*inv_Kii_l{s}*Kib_l{s};
            bp_l{s} = fb_l{s}-Kbi_l{s}*inv_Kii_l{s}*fi_l{s};
            t_init_local = toc(t_start);
        else
    
            Kii_l{s} = K_l{s}(i_list{s}, i_list{s});
            Kib_l{s} = K_l{s}(i_list{s}, b_list{s});
            Kbi_l{s} = K_l{s}(b_list{s}, i_list{s});
            Kbb_l{s} = K_l{s}(b_list{s}, b_list{s});
        
            fi_l{s} = f_l{s}(i_list{s});
            fb_l{s} = f_l{s}(b_list{s});
    
            inv_Kii_l{s} = inv(Kii_l{s});
            Sp_l{s} = Kbb_l{s}-Kbi_l{s}*inv_Kii_l{s}*Kib_l{s};
            bp_l{s} = fb_l{s}-Kbi_l{s}*inv_Kii_l{s}*fi_l{s};
        end
    
        
    end
    Sp_diam = blkdiag(Sp_l{:}); % diagonal concatenation
    bp_diam = cat(1,bp_l{:});   % 1 for concatenation in first dim : vertical    
    A_diam = A_op(N);
    Ab_diam = A_bar_op(N);
    Sp = A_diam*Sp_diam*A_diam';
    bp = A_diam*bp_diam;

    % solution of primal shur complement with direct method
    t_start = tic;
    ub = Sp\bp;
    t_sp = toc(t_start);
    % ---------------------------------------------------------------------
    %% Distributed conjugate gradient to solve primal interface problem----
        
    %Initialization phase
    m = N;
    epsilon = 1e-10;
    ub = zeros(N-1, 1); %% n rows = n of boundary nodes
    ub_diam = Ab_diam'*ub;
    
    ub_l=vert_diam_to_s(ub_diam, N);
    
    %solve the dirichlet problem 
    ui_l=cell(1,N);
    
    rb_l=cell(1,N);
    for s=1:N
            if s==1
                t_start_cg = tic;
            end
            ui_l{s} = inv_Kii_l{s} * (fi_l{s} - Kib_l{s}*ub_l{s});
            
            Kbi_bb = cat(2, Kbi_l{s}, Kbb_l{s});   
    
            uib = cat(1, ui_l{s}, ub_l{s});
            
            rb_l{s} = -( Kbi_bb * uib - fb_l{s});

            if s==1
                t_init_cg = toc(t_start_cg);
            end
        
    end
     %time_taken_dist_primal_CG(j) = t_init_cg;
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
            if s==1
                t_start_cg = tic;
            end
            di_l{s} = -inv_Kii_l{s}*Kib_l{s}*db_l{s};
            
            %compute local matrix vector product
            Kbi_bb = cat(2, Kbi_l{s}, Kbb_l{s});
            dib = cat(1, di_l{s}, db_l{s});
            
            Sp_db_l{s} = Kbi_bb*dib;
            s;
            Sp_db_l{s};
            if s==1
                t_local_cg = toc(t_start_cg);
            end
            
        end
        
        % Conjugate gradient update step
        
        t_loop = tic;
        
        Sp_db_diam = cat(1, Sp_db_l{:});
        Sp_db = A_diam*Sp_db_diam;             %Compute global matrix-vector product
        
        alpha = (rb'*rb) / (db'*Sp_db);        %Compute optimal step
       
        ub_next = ub + alpha*db;               %Compute iterate
    
        rb_next = rb - alpha*Sp_db;            %Compute residual
        
        beta =  (rb_next'*rb_next) / (rb'*rb); %Compute orthogonalization parameter 
        
        db = rb_next + beta*db;                %Update search direction
    
        rb = rb_next;
        ub = ub_next;
        
        t_loop_cg = toc(t_loop);

        time_taken_dist_primal_CG(j) =time_taken_dist_primal_CG(j) + t_loop_cg;
        %%%--------convergence criterion-----
        if norm(rb)/r0_norm < epsilon 
            break;
        end
        %%%----------------------------------
        

    end
    
    ub_diam=A_diam'*ub;
    ub_l=vert_diam_to_s(ub_diam, N);
    
    ui_l = cell(1,N);
    for s=1:N
        if s == 1        
            Kii = K_l{s}(i_list{s}(2:end), i_list{s}(2:end));
            Kib = K_l{s}(i_list{s}(2:end), b_list{s});
            fi = f_l{s}(i_list{s}(2:end));
    
            t_start = tic;
            ui_l{s} = Kii\fi-Kib*ub_l{s};
            t_ui = toc(t_start);
    
        else
            Kii = K_l{s}(i_list{s}, i_list{s});
            Kib = K_l{s}(i_list{s}, b_list{s});
            fi = f_l{s}(i_list{s});
    
            
            ui_l{s} = inv(Kii)*(fi-Kib*ub_l{s});

        end
        %ui_l{s};
    end
    time_taken_primal(j) = t_sp;
    time_taken_dist_primal_CG(j) =time_taken_dist_primal_CG(j);
    j=j+1;
end

%semilogy(n_min_dofs:n_min_dofs:n_max_dofs, time_taken_fem, Color="blue")
%hold on
semilogy(n_min_dofs:n_min_dofs:n_max_dofs, time_taken_primal, Color="green")
hold on
semilogy(n_min_dofs:n_min_dofs:n_max_dofs, time_taken_dist_primal_CG, Color="red")

legend("Primal direct invert Sp", "Distributed CG")

xlabel("Total number of drees of freedom")
ylabel("Time in logaritmic scale log s")
title("Runtime comparison ")
hold off


