clear all; clc; close all;

% --- Paramètres fixes ---
L = 30; 
S = 1; 
E = 2e5; 
Fd = 10;
tol = 1e-6; 

% On fixe le nombre d'éléments par sous-domaine pour éviter les problèmes d'arrondi
% Cela assure que le maillage est toujours valide quel que soit H
n_elements_par_SD = 10; 

% Plage de N (Nombre de sous-domaines) à tester
N_values = 2:30; 

results_iter = [];

fprintf('--- Étude de la N-Scalabilité (Variation du nombre de sous-domaines) ---\n');
fprintf('%-10s %-10s %-15s\n', 'N', 'H', 'Iterations');

for N = N_values
    H = L / N;              % Taille du sous-domaine
    h = H / n_elements_par_SD; % Taille de l'élément (adaptée pour tomber juste)
    
    % Appel de la fonction (Utilise ta fonction corrigée primal_CJ_f)
    [~, k] = primal_CJ_f(L, S, E, Fd, H, h, tol);
    
    results_iter = [results_iter, k];
    
    fprintf('%-10d %-10.4f %-15d\n', N, H, k);
end

% --- Tracé du graphique ---
figure('Name', 'N-Scalabilité Primal Schur', 'Color', 'w');

plot(N_values, results_iter, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'Color', 'r');

xlabel('Nombre de sous-domaines (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Nombre d''itérations pour converger', 'FontSize', 12, 'FontWeight', 'bold');
title({'Scalabilité par rapport au nombre de sous-domaines', '(Primal Schur sans préconditionneur)'}, 'FontSize', 14);
grid on;

% Ajout d'une annotation pour l'interprétation
%text(mean(N_values), mean(results_iter), 'Augmentation \rightarrow Non Scalable', ...
%    'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w');