function [i_list,b_list] = node_lists(N, n)
%NODE_LISTS generates 2 lists: internal nodes and boundary nodes lists
%   i_list internal nodes numbering for each substructure
%   b_list interface nodes numbering for each substructure
arguments (Input)
    N
    n
end

arguments (Output)
    i_list
    b_list
end

i_list = cell(1,N); %internal nodes 
b_list = cell(1,N); %interface nodes

for s=1:N % for each substructure
    for e=1:n+1          % for each node in substructure
        if e == 1
            if s == 1 % first node of first substructure is internal 
                i_list{s} = [i_list{s} e];
            else
                b_list{s} = [b_list{s} e]; %otherwise it's interface
            end
        elseif e == n+1
            if s == N % last node of first substructure is internal 
                i_list{s} = [i_list{s} e];
            else
                b_list{s} = [b_list{s} e]; %otherwise it's interface
            end
        else
            i_list{s} = [i_list{s} e]; %otherwise it's an internal node
        end
    end 
end

end