function [K_l,f_l] = global_k_f_lists(N, H, E, S, Fd, n, node)
%GLOBAL_K_F_LISTS Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    N
    H
    E
    S
    Fd
    n
    node
end

arguments (Output)
    K_l
    f_l
end

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
end