clear all; clc; close all;

% --- Paramètres du problème ---
L = 30; 
S = 1; 
E = 2e5; 
Fd = 10; % (Note: Fd n'influence pas le conditionnement, mais on le garde pour la forme)
H = 10;  % Taille fixe des sous-domaines (Donc N = 3 sous-domaines)

% --- Paramètres de l'étude (Question 2.2) ---
% On fait varier n (nombre d'éléments par sous-domaine)
% Comme h = H/n, cela fait varier h.
n_values = [2, 4, 8, 16, 32, 60]; 

% Tableaux pour stocker les résultats
store_h_H = [];   % Rapport h/H
store_h = [];     % Taille h
store_cond_Sp = []; % Conditionnement de Sp
store_cond_K = [];  % Conditionnement de K global

disp('Début de l''analyse de scalabilité...');
fprintf('%-10s %-10s %-15s %-15s\n', 'n', 'h', 'cond(Sp)', 'cond(K)');

for n = n_values
    h = H / n; % Taille de l'élément
    
    % Appel de la fonction séparée (Fichier 1)
    [c_Sp, c_K] = calcul_cond_2_2(L, S, E, H, h);
    
    % Stockage
    store_h_H = [store_h_H, h/H];
    store_h = [store_h, h];
    store_cond_Sp = [store_cond_Sp, c_Sp];
    store_cond_K = [store_cond_K, c_K];
    
    % Affichage dans la console
    fprintf('%-10d %-10.4f %-15.4e %-15.4e\n', n, h, c_Sp, c_K);
end

% --- Tracé des courbes ---

% Figure pour Sp (Question 2.2a)
figure(1);
loglog(store_h_H, store_cond_Sp, '-o', 'LineWidth', 2);
xlabel('h / H (finesse relative)', 'FontSize', 12);
ylabel('Conditionnement \kappa(S_p)', 'FontSize', 12);
title('Conditionnement de la matrice de Schur Primal', 'FontSize', 14);
grid on;

% Figure pour K global (Question 2.2b)
figure(2);
loglog(store_h, store_cond_K, '-r+', 'LineWidth', 2);
xlabel('h (taille élément)', 'FontSize', 12);
ylabel('Conditionnement \kappa(K)', 'FontSize', 12);
title('Conditionnement de la matrice de rigidité Globale', 'FontSize', 14);
grid on;

disp('Analyse terminée. Figures générées.');