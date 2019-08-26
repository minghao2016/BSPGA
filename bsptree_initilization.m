function [bsptree, cur_tree_size] = bsptree_initilization(dim, max_tree_size)
 bsptree = struct (       'x',              zeros(max_tree_size,  dim),                     ...% SOLUTION %
                          'optimal',         zeros(1,dim),                                   ...% OPTIMAL SOLUTION
                          'min_f',          Inf,                                            ...% MIN FITNESS VALUE %
                          'cro_d',          -1,                                             ...% CROSSOVER DIM %
                          'f',              zeros(max_tree_size,  1  ),                     ...% FITNESS OF SOLUTION %
                          'd',              -1 * ones(max_tree_size,1),                     ...% PARTITIONING ON THIS DIMENSION %
                          'p_idx',          zeros(max_tree_size,   1 ),                     ...% POINTER TO PARENT %
                          'l_idx',          (max_tree_size + 1) * ones(max_tree_size, 1),	...% POINTER TO LEFT CHILD %
                          'r_idx',          (max_tree_size + 1) * ones(max_tree_size, 1)	...% POINTER TO RIGHT CHILD %
                          );
            bsptree.d(1, 1) = -2;
            cur_tree_size = 0;
end                                                                                                                                                                                                                                                                                                                                                             