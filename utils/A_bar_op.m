function [Ab_diam] = A_bar_op(N)
%A_BAR_OP Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    N
end

arguments (Output)
    Ab_diam  
end

% build A bar diamond
Ab_l = cell(1,N);
for s=1:N
    if s == 1 || s == N
        Ab_l{s} = zeros(N-1, 1); 
    else
        Ab_l{s} = zeros(N-1, 2);
    end
end

for c=1:(N-1)
    if c == 1
        Ab_l{c}(c, 1) = 1;
        Ab_l{c+1}(c, 1) = -1; 
    else
        Ab_l{c}(c, 2) = 1; 
        Ab_l{c+1}(c, 1) = -1; 
    end
    
end

Ab_diam = cat(2, Ab_l{:}); % 2: horizontal concatenation
end