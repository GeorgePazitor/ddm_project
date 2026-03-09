%% --- POST-TRAITEMENT : Question 2.9 ---

ub_l = vert_diam_to_s(A_diam'*ub, N); % On les remet sous forme de liste par sous-domaine

% 4. Calcul des noeuds internes et assemblage pour le tracé
U_total = zeros(N*n + 1, 1);
coord_x = linspace(0, L, N*n + 1)';

for s = 1:N
    % Calcul des déplacements internes ui_s via le problème de Dirichlet local
    ui_s = Kii_l{s} \ (fi_l{s} - Kib_l{s} * ub_l{s});
    
    global_nodes = (s-1)*n + (1 : n+1);
    
    % On place les interfaces
    U_total(global_nodes(b_list{s})) = ub_l{s};
    
    % On place les noeuds internes (avec gestion de l'encastrement)
    if length(ui_s) == length(i_list{s}) - 1
        % Si ui_s est plus court, c'est que le noeud 1 (encastré) a été retiré de la matrice !
        U_total(global_nodes(i_list{s}(2:end))) = ui_s;
        U_total(global_nodes(i_list{s}(1))) = 0; % On force le 0 de l'encastrement
    else
        % Cas normal (sous-domaines flottants)
        U_total(global_nodes(i_list{s})) = ui_s;
    end
end

% 5. Plot (Requested plot)
figure('Name', 'Solution FETI 2.9', 'Color', 'black');
plot(coord_x, U_total, 'r.', 'MarkerSize', 10, 'DisplayName', 'DDM Solution (FETI)');
hold on;
plot(coord_x, (Fd/(E*S))*coord_x, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Exact Analytical Solution', 'Color', 'green');
xlabel('Position x (mm)', 'FontSize', 12, 'Color', 'w');
ylabel('Displacement u(x) (mm)', 'FontSize', 12, 'Color', 'w');
legend('Location', 'northwest', 'FontSize', 11, 'Color', 'black');
title('Reconstruction of Displacements after FETI', 'FontSize', 14, 'Color', 'w');
grid on;