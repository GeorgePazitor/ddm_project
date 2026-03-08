clear all; clc; close all;

% --- Paramètres fixes ---
L = 30; 
S = 1; 
E = 2e5; 
Fd = 10;
tol = 1e-6; 

% On fixe le nombre d'éléments par sous-domaine
n_elements_par_SD = 10; 

% Plage de N (Nombre de sous-domaines) à tester
N_values = 2:30; 

results_iter = [];

fprintf('--- Étude de la N-Scalabilité BDD (Variation du nombre de sous-domaines) ---\n');
fprintf('%-10s %-10s %-15s\n', 'N', 'H', 'Iterations');

for N = N_values
    H = L / N;              % Taille du sous-domaine
    h = H / n_elements_par_SD; 
    
    % Appel de la fonction BDD (Préconditionnée)
    [~, k] = primal_precond_CJ_f(L, S, E, Fd, H, h, tol);
    
    results_iter = [results_iter, k];
    
    fprintf('%-10d %-10.4f %-15d\n', N, H, k);
end

% --- Plot graph ---
figure('Name', 'N-Scalability BDD', 'Color', 'black');

plot(N_values, results_iter, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'g', 'Color', 'g');

xlabel('Number of subdomains (N)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w');
ylabel('Number of iterations to converge', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w');
title({'Scalability vs Number of subdomains N', '(BDD method: Neumann-Neumann + Coarse)'}, 'FontSize', 14, 'Color', 'w');
grid on;
ylim([0 max(results_iter)+5]);

% Annotation
text(mean(N_values), mean(results_iter)+0.5, 'Plateau \rightarrow SCALABLE!', ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w');