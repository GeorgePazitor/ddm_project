function [A_diam] = A_op(N)
%A_OP build the primal assembly operator A
%   Based on the knowledge of the geometry of the 1D problem 
arguments (Input)
    N
end

arguments (Output)
    A_diam
end

% build A diamond
A_l = cell(1,N);        
for s=1:N
    %for the first and the last substruct it's just a vector A(s) 
    % has 1 column (1 single interaction per substructure)
    if s == 1 || s == N   
        A_l{s} = zeros(N-1, 1); 
    %for the first and the last substruct it's just a vector A(s) has 2 col
    else                 
        A_l{s} = zeros(N-1, 2);
    end
end

for c=1:(N-1)
    if c == 1
        A_l{c}(c, 1) = 1;
        A_l{c+1}(c, 1) = 1; 
    else
        A_l{c}(c, 2) = 1; 
        A_l{c+1}(c, 1) = 1; 
    end
end

A_diam = cat(2, A_l{:}); 
end