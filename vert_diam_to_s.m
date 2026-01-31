function [v_l] = vert_diam_to_s(v_diam, N)
%VERT_DIAM_TO_S converts a convatenated version of a vector to a list in s
%   The concatenation of substructure specific vectors are identified by 
%   the diamond, this function take the single vectors from the concatenat.
%   and place it in the corresponding position of the substructure.
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
