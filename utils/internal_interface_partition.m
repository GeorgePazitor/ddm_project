function [Kii_l, Kib_l, Kbi_l, Kbb_l, fi_l, fb_l] = internal_interface_partition(K_l, f_l, i_list, b_list, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    K_l
    f_l
    i_list
    b_list
    N
end

arguments (Output)
    Kii_l
    Kib_l
    Kbi_l
    Kbb_l
    fi_l
    fb_l
end

Kii_l = cell(1, N);
Kib_l = cell(1, N);
Kbi_l = cell(1, N);
Kbb_l = cell(1, N);

fi_l = cell(1, N);
fb_l = cell(1, N);

for s=1:N    
    if s == 1 %elimination of the degree of freedom to impose ud
        %just exclude the first line and/or column for internal dofs 
        Kii_l{s} = K_l{s}(i_list{s}(2:end), i_list{s}(2:end));
        Kib_l{s}= K_l{s}(i_list{s}(2:end), b_list{s});
        Kbi_l{s} = K_l{s}(b_list{s}, i_list{s}(2:end));
        Kbb_l{s} = K_l{s}(b_list{s}, b_list{s});

        fi_l{s} = f_l{s}(i_list{s}(2:end));
        fb_l{s} = f_l{s}(b_list{s});
    else

        Kii_l{s} = K_l{s}(i_list{s}, i_list{s});
        Kib_l{s} = K_l{s}(i_list{s}, b_list{s});
        Kbi_l{s} = K_l{s}(b_list{s}, i_list{s});
        Kbb_l{s} = K_l{s}(b_list{s}, b_list{s});
    
        fi_l{s} = f_l{s}(i_list{s});
        fb_l{s} = f_l{s}(b_list{s});
    end


end
end