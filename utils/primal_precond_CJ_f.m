function [U_total, n_iter] = primal_precond_CJ_f(L, S, E, Fd, H, h, tol)
    % PRIMAL_PRECOND_CJ_F : Méthode BDD (Neumann-Neumann + Problème Grossier)
    
    %% 1. Initialisation
    if nargin < 1 || isempty(H), H = 10; end
    if nargin < 2 || isempty(h), h = 5; end
    if nargin < 3 || isempty(L), L = 30; end
    if nargin < 4 || isempty(S), S = 1; end
    if nargin < 5 || isempty(Fd),Fd = 10; end
    if nargin < 6 || isempty(E), E = 2e5; end
    if nargin < 7 || isempty(tol), tol = 1e-6; end
    n = round(H/h); 
    N = round(L/H); 
    k0 = E * S / h; 

    A_diam = A_op(N); % Opérateur d'assemblage
    
    %% 2. Construction des Matrices Locales et Modes Rigides
    % On n'utilise plus node_lists standard car on doit gérer finement les tailles
    
    K_l = cell(1,N); 
    f_l = cell(1,N);
    Rb_l = cell(1,N);
    
    % Variables pour le Schur
    Sp_l = cell(1,N);
    bp_l = cell(1,N);
    
    % Pour préconditionneur
    inv_Kii_l = cell(1,N);

    for s=1:N 
        % --- A. Matrice de rigidité locale complète ---
        K_local = k0 * (diag(2*ones(n+1,1)) - diag(ones(n,1),1) - diag(ones(n,1),-1));
        K_local(1,1) = k0; K_local(n+1,n+1) = k0;
        
        f_local = zeros(n+1,1);
        f_local(end) = Fd; % Force à droite (simplifié)
        
        % --- B. Gestion des Interfaces et Modes Rigides ---
        if s == 1
            % Fixé à gauche (Dirichlet)
            % Interface : Droite seulement (noeud n+1)
            K_local(1,:) = 0; K_local(:,1) = 0; K_local(1,1) = 1; f_local(1) = 0;
            
            b_list = n+1;       % Interface (1 noeud)
            i_list = 1:n;       % Internes (inclut le noeud 1 bloqué)
            
            % Pas de mode rigide car fixé
            Rb_l{s} = zeros(1, 0); 
            
        elseif s == N
            % Libre à droite
            % Interface : Gauche seulement (noeud 1)
            % Le noeud n+1 est "interne" au sens de Schur (il est condensé)
            b_list = 1;         % Interface (1 noeud)
            i_list = 2:n+1;     % Internes
            
            % Mode rigide : Translation [1]. (1x1)
            Rb_l{s} = 1; 
            
        else
            % Flottant au milieu
            % Interfaces : Gauche (1) et Droite (n+1)
            b_list = [1, n+1];  % Interface (2 noeuds)
            i_list = 2:n;       % Internes
            
            % Mode rigide : Translation [1; 1] (normalisée)
            Rb_l{s} = [1; 1] / sqrt(2);
        end
        
        % Stockage K et f
        K_l{s} = K_local;
        f_l{s} = f_local;
        
        % --- C. Calcul du Schur local Sp ---
        K_ii = K_local(i_list, i_list);
        K_bb = K_local(b_list, b_list);
        K_ib = K_local(i_list, b_list);
        K_bi = K_local(b_list, i_list);
        
        inv_Kii_l{s} = inv(K_ii);
        
        Sp_l{s} = K_bb - K_bi * (inv_Kii_l{s} * K_ib);
        
        % Second membre condensé
        fi = f_local(i_list);
        fb = f_local(b_list);
        bp_l{s} = fb - K_bi * (inv_Kii_l{s} * fi);
    end

    %% 3. Assemblage des Opérateurs Globaux
    Sp_diam = blkdiag(Sp_l{:});
    Sp = A_diam * Sp_diam * A_diam';
    
    bp_diam = cat(1, bp_l{:});
    bp = A_diam * bp_diam;
    
    Rb_diam = blkdiag(Rb_l{:}); % Doit maintenant avoir la bonne taille !

    %% 4. Construction du Problème Grossier (Coarse Problem)
    % Matrice de pondération D (Partition de l'unité)
    % Pour faire simple : D_inv = diag(multiplicité)
    mult = A_diam * ones(size(A_diam,2), 1); 
    D_inv = spdiags(1./mult, 0, length(mult), length(mult));
    
    % A_scaled correspond à "A * D" (moyenne aux interfaces)
    A_scaled = D_inv * A_diam; 
    
    % G : Base des modes rigides assemblés
    G = A_scaled * Rb_diam; 
    
    % Opérateur du problème grossier S0
    S0 = G' * Sp * G;
    
    %% 5. Algorithme BDD (Gradient Conjugué Préconditionné + Projeté)
    
    % Solution initiale (Problème grossier)
    lambda_0 = S0 \ (G' * bp);
    u = G * lambda_0; 
    
    r = bp - Sp * u; % Résidu initial
    
    % Fonction locale Neumann
    solve_neumann = @(res) apply_neumann(res, N, A_scaled, Sp_l);
    
    % Préconditionnement initial
    z = solve_neumann(r);
    
    % Projection initiale (d = P*z)
    % P*z = z - G * S0^-1 * G' * Sp * z
    correction = S0 \ (G' * Sp * z);
    d = z - G * correction;
    
    rz = r' * d; % Produit scalaire <r, d> (d est déjà precond+projeté)
    
    n_iter = 0;
    max_iter = 200;
    
    for k = 1:max_iter
        % Produit matrice-vecteur
        w = Sp * d; 
        
        % Pas de descente alpha
        if (d' * w) == 0, break; end % Sécurité
        alpha = rz / (d' * w);
        
        % Mise à jour
        u = u + alpha * d;
        r = r - alpha * w;
        
        if norm(r) < tol * norm(bp) % Critère relatif souvent mieux
            n_iter = k;
            break;
        end
        if norm(r) < tol % Critère absolu si bp petit
            n_iter = k;
            break;
        end
        
        % Préconditionnement
        z_new = solve_neumann(r);
        
        % Projection de la direction (Problème grossier)
        correction_new = S0 \ (G' * Sp * z_new);
        z_projected = z_new - G * correction_new;
        
        % Beta
        rz_new = r' * z_projected; 
        beta = rz_new / rz;
        
        % Direction suivante
        d = z_projected + beta * d;
        
        rz = rz_new;
        n_iter = k; % Mise à jour itération courante
    end
    
    U_total = zeros(n*N+1, 1); % Placeholder
end

function z_global = apply_neumann(r_global, N, A_scaled, Sp_l)
    % 1. Distribuer le résidu (r_loc = D * A' * r_glob = A_scaled' * r_glob)
    r_local_all = A_scaled' * r_global;
    
    % 2. Résoudre localement (Neumann)
    z_parts = cell(1, N);
    current_idx = 1;
    
    for s = 1:N
        n_dof_interface = size(Sp_l{s}, 1);
        idx_range = current_idx : (current_idx + n_dof_interface - 1);
        
        r_s = r_local_all(idx_range);
        
        % Résolution pseudo-inverse (car modes rigides possibles)
        z_parts{s} = pinv(full(Sp_l{s})) * r_s;
        
        current_idx = current_idx + n_dof_interface;
    end
    
    % 3. Assembler (z_glob = A * D * z_loc = A_scaled * z_loc)
    z_diam = cat(1, z_parts{:});
    z_global = A_scaled * z_diam;
end