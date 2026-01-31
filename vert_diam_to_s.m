function [v_l] = vert_diam_to_s(v_diam, N)
%VERT_DIAM_TO_S Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    v_diam
    N    
end

arguments (Output)
    v_l
end

v_l=cell(1,N);

for s=1:N
    if s == 1
        v_l{s} = v_diam(s);
    else
        if s==N
            v_l{s} = v_diam((s*2)-2);
        else
            v_l{s} = v_diam(((s*2)-2):((s*2)-1));
        end
    end
end
