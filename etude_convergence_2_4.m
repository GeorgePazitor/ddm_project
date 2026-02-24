clear all; clc; close all;

% --- Paramètres fixes ---
L = 30; 
S = 1; 
E = 2e5; 
Fd = 10;
tol = 1e-6; 

% Liste des finesses (n = H/h) à tester
n_values = 2.^(1:8); %[2, 4, 8, 16, 32, 64, 128]; 

%% ÉTUDE 1 : Cas de l'exercice (H = 10 => N = 3 sous-domaines)
H1 = 10;
results_h_H1 = [];
results_iter1 = [];

fprintf('--- Début Étude 1 (H=%d, N=%d) ---\n', H1, L/H1);
fprintf('%-10s %-10s %-15s\n', 'n', 'h/H', 'Iterations');

for n = n_values
    h = H1 / n;
    
    % Appel de ta fonction corrigée
    [~, k] = primal_CJ_f(L, S, E, Fd, H1, h, tol);
    
    results_h_H1 = [results_h_H1, h/H1];
    results_iter1 = [results_iter1, k];
    
    fprintf('%-10d %-10.4f %-15d\n', n, h/H1, k);
end

%% ÉTUDE 2 : Cas "Démonstration" (H = 1 => N = 30 sous-domaines)
H2 = 1; 
results_h_H2 = [];
results_iter2 = [];

fprintf('\n--- Début Étude 2 (H=%d, N=%d) ---\n', H2, L/H2);
fprintf('%-10s %-10s %-15s\n', 'n', 'h/H', 'Iterations');

for n = n_values
    h = H2 / n;
    
    % Appel de ta fonction corrigée
    [~, k] = primal_CJ_f(L, S, E, Fd, H2, h, tol);
    
    results_h_H2 = [results_h_H2, h/H2];
    results_iter2 = [results_iter2, k];
    
    fprintf('%-10d %-10.4f %-15d\n', n, h/H2, k);
end

%% Tracé des graphiques

figure('Name', 'Etude de Scalabilité - Primal Schur', 'Color', 'w');

% Sous-graphe 1 : Cas de l'exercice
subplot(1, 2, 1);
semilogx(results_h_H1, results_iter1, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('h / H (Finesse relative)');
ylabel('Nombre d''itérations');
title(['Cas Exercice : N = ', num2str(L/H1), ' sous-domaines']);
grid on;
ylim([0 max(results_iter1)+5]); % Ajuste l'échelle
set(gca, 'XDir', 'reverse'); % On inverse l'axe X pour voir h/H diminuer vers la droite (raffinement)

% Sous-graphe 2 : Cas Démonstration
subplot(1, 2, 2);
semilogx(results_h_H2, results_iter2, '-ro', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('h / H (Finesse relative)');
ylabel('Nombre d''itérations');
title(['Cas Démonstration : N = ', num2str(L/H2), ' sous-domaines']);
grid on;
set(gca, 'XDir', 'reverse');

sgtitle('Question 2.4 : Scalabilité du Gradient Conjugué (Sans Précond.)');