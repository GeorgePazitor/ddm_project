function [cond_Sp, cond_K_global] = calcul_cond_2_2(L, S, E, H, h)
    % CALCUL_COND Calcule le conditionnement de Sp et de K global
    
    %% initialisation of parameters
    if nargin < 1 || isempty(H), H = 10; end
    if nargin < 2 || isempty(h), h = 5; end
    if nargin < 3 || isempty(L), L = 30; end
    if nargin < 4 || isempty(S), S = 1; end
    %if nargin < 5 || isempty(Fd),Fd = 10; end
    if nargin < 6 || isempty(E), E = 2e5; end
    %if nargin < 7 || isempty(tol), tol = 1e-6; end

    %% 1. Paramètres dérivés
    N = round(L/H);     % Nombre de sous-domaines
    n = round(H/h);     % Nombre d'éléments par sous-domaine
    k0 = E * S / h;     % Raideur élémentaire
    
    %% 2. Calcul du Complément de Schur Primal (Sp)
    Sp_l = cell(1, N); 

    % Boucle sur chaque sous-domaine
    for s = 1:N
        % Construction K_l
        K_l = sparse(n+1, n+1);
        K_l = k0 * (diag(2*ones(n+1,1)) - diag(ones(n,1),1) - diag(ones(n,1),-1));
        K_l(1,1) = k0; K_l(n+1,n+1) = k0;
        
        % --- CORRECTION ICI : Gestion précise des noeuds ---
        if s == 1
            % Sous-domaine 1 : 
            % - Gauche (noeud 1) : Bloqué (Dirichlet), on l'enlève.
            % - Droite (noeud n+1) : Interface.
            i_nodes = 2:n;      % Noeuds internes
            b_nodes = n+1;      % Noeud d'interface
            
        elseif s == N
            % Sous-domaine N (Dernier) :
            % - Gauche (noeud 1) : Interface.
            % - Droite (noeud n+1) : Libre (Effort imposé), donc c'est un noeud INTERNE pour Schur.
            i_nodes = 2:(n+1);  % On inclut le dernier noeud dans les internes à condenser
            b_nodes = 1;        % Seul le noeud de gauche est sur l'interface
            
        else
            % Sous-domaines intermédiaires :
            % - Gauche (1) et Droite (n+1) sont des interfaces.
            i_nodes = 2:n;
            b_nodes = [1, n+1];
        end
        % ---------------------------------------------------
        
        % Schur local
        K_ii = K_l(i_nodes, i_nodes);
        K_bb = K_l(b_nodes, b_nodes);
        K_ib = K_l(i_nodes, b_nodes);
        K_bi = K_l(b_nodes, i_nodes);
        
        Sp_l{s} = K_bb - K_bi * (K_ii \ K_ib);
    end
    
    % Assemblage global de Sp
    A = A_op(N);                
    Sp_diam = blkdiag(Sp_l{:}); 

    % Calcul de Sp assemblé
    Sp = A * Sp_diam * A';      
    
    % Calcul du conditionnement de Sp
    cond_Sp = cond(full(Sp));
    
    %% 3. Calcul de la rigidité globale K
    n_total_elem = N * n;
    n_total_nodes = n_total_elem + 1;
    
    K_global = sparse(n_total_nodes, n_total_nodes);
    K_global = k0 * (diag(2*ones(n_total_nodes,1)) - diag(ones(n_total_elem,1),1) - diag(ones(n_total_elem,1),-1));
    K_global(1,1) = k0; 
    K_global(end,end) = k0;
    
    % Condition limite u(0)=0 : on enlève la première ligne/colonne
    K_reduit = K_global(2:end, 2:end);
    
    cond_K_global = cond(full(K_reduit));
end