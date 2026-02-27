clear all; clc; close all;

% --- Paramètres fixes ---
L = 30; S = 1; E = 2e5; Fd = 10;
H = 10; % Taille de sous-domaine fixe (N=3)

% Liste des finesses à tester
n_values = [2, 4, 8, 16, 32, 64, 128]; 

results_h_H = [];
results_cond_Prec = [];

fprintf('--- Question 2.6 : Conditionnement du système préconditionné ---\n');
fprintf('%-10s %-10s %-15s\n', 'n', 'h/H', 'cond(M^-1 * Sp)');

for n = n_values
    h = H / n;
    
    % --- 1. Construction des matrices (inspiré de primal_precond_CJ_f) ---
    N = round(L/H); 
    k0 = E * S / h; 
    A_diam = A_op(N); 
    
    % Matrices locales
    Sp_l = cell(1,N);
    for s=1:N 
        % Matrice K locale
        K_local = k0 * (diag(2*ones(n+1,1)) - diag(ones(n,1),1) - diag(ones(n,1),-1));
        K_local(1,1) = k0; K_local(n+1,n+1) = k0;
        
        % Gestion des interfaces
        if s == 1
            i_list = 1:n; b_list = n+1;
            K_local(1,:) = 0; K_local(:,1) = 0; K_local(1,1) = 1; % Dirichlet
        elseif s == N
            i_list = 2:n+1; b_list = 1;
        else
            i_list = 2:n; b_list = [1, n+1];
        end
        
        K_ii = K_local(i_list, i_list);
        K_bb = K_local(b_list, b_list);
        K_ib = K_local(i_list, b_list);
        K_bi = K_local(b_list, i_list);
        
        % Schur local
        Sp_l{s} = K_bb - K_bi * (K_ii \ K_ib);
    end
    
    % Assemblage de Sp
    Sp_diam = blkdiag(Sp_l{:});
    Sp = A_diam * Sp_diam * A_diam';
    
    % --- 2. Construction du Préconditionneur Neumann-Neumann (M^-1) ---
    % D_inv : Matrice de pondération (Partition de l'unité)
    mult = A_diam * ones(size(A_diam,2), 1); 
    D_inv = spdiags(1./mult, 0, length(mult), length(mult));
    A_scaled = D_inv * A_diam; 
    
    % Inverses locaux (Pseudo-inverses car modes rigides possibles)
    inv_Sp_l = cell(1,N);
    for s=1:N
        inv_Sp_l{s} = pinv(full(Sp_l{s}));
    end
    
    % Assemblage du préconditionneur : M^-1 = A_scaled * blockdiag(Sp_inv) * A_scaled'
    inv_Sp_diam = blkdiag(inv_Sp_l{:});
    M_inv = A_scaled * inv_Sp_diam * A_scaled';
    
    % --- 3. Calcul du conditionnement du produit ---
    % Matrice préconditionnée
    P_mat = M_inv * Sp;
    
    % Conditionnement
    c = cond(full(P_mat));
    
    results_h_H = [results_h_H, h/H];
    results_cond_Prec = [results_cond_Prec, c];
    
    fprintf('%-10d %-10.4f %-15.4f\n', n, h/H, c);
end

% --- Tracé ---
figure('Name', 'Conditionnement Préconditionné', 'Color', 'black');
semilogx(results_h_H, results_cond_Prec, '-go', 'LineWidth', 2);
xlabel('h / H (Finesse relative)', 'Color', 'w');
ylabel('Conditionnement k(tilde{S}_p^{-1} S_p)', 'Color', 'w');
title('Question 2.6 : Scalabilité du Conditionnement brut (O(1))');
grid on;
set(gca, 'XDir', 'reverse'); % Axe inversé pour montrer le raffinement vers la droite
ylim([0 max(results_cond_Prec)*1.2]);