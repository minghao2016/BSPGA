% ===   ***** ***** ***** ***** ***** ***** *****  	=== %
% ===   BSPTREE_INSERTION                         	=== %
% ===   ***** ***** ***** ***** ***** ***** *****  	=== %
% === 	CREATED: 22 NOVEMBER, 2014                 	=== %
% ===   UPDATED: 23 MAY, 2015                       === %
% ===   ***** ***** ***** ***** ***** ***** *****  	=== %

function [cur_bsptree_node, x, cur_popf, cro_d] = bsptree_insertion(x, cur_bsptree_node,gen)
global BSPTREE;
global test;
testnumber = 1;                  % TEST FUCTION'S NUMBER
umin = -10;                      % LOWER BOUND
umax = 10;                       % UPPER BOUND
D = 40;                          % FUNCTION DIMENSION
Length = 12;                     % xx LENGTH
[pop_size, dim] = size(x);
cur_popf = Inf * ones(pop_size, 1);
for idx = 1 : pop_size
    if BSPTREE.d(1, 1) == -2  % ROOT NODE ONLY %
        BSPTREE.x(1,: ) = x(1,:);
        
        %--------------UNCODING------------------%
        y = zeros(D,1);
        m = inf*ones(D,Length);
        for i = 1:dim
            if mod(i,Length) == 0
                m(i/Length,:) = x(1,i-Length+1:i);
            end
        end
        for i = 1:D
            for j = 1:1:Length
                y(i) = y(i) + m(i,j) * 2^(j-1);
            end
            xx(i) = (umax-umin)*y(i)/4095+umin;
        end
        %-----------END UNCODING-----------------%
        
        BSPTREE.f(1,:) =  Test_Function(xx, D, testnumber);
        
        cur_popf(1, 1) = BSPTREE.f(1, 1);
        BSPTREE.optimal = BSPTREE.x(1,: );
        BSPTREE.min_f = BSPTREE.f(1,1);
        BSPTREE.d(1, 1) = 0;   % THE ROOT IS NOW A 'LEAF' %
        cur_bsptree_node = 1;
    else  % MEANS IF 'BSPTREE.D(1, 1) ~= -2' %
        n_idx = 1;
        flag = 0;
        h = 1;
        testCount = 0;
        while BSPTREE.d(n_idx, 1)
            testCount = testCount+1;
            jdx = BSPTREE.d(n_idx, 1);
            if x(idx,jdx) > 0.5
                n_idx = BSPTREE.l_idx(n_idx, 1);
                memo_x = BSPTREE.x(n_idx,:);
                if all(memo_x == x(idx,:)) && flag == 0
                    flag = 1;
                end
            else
                n_idx = BSPTREE.r_idx(n_idx, 1);
                memo_x = BSPTREE.x(n_idx,:);
                if all(memo_x == x(idx,:)) && flag == 0
                    flag = 1;
                end
            end
            h = h + 1;
            if h > dim
                h = h - 1;
                break;
            end
        end
        test(gen,idx) = testCount;
        memo_x = BSPTREE.x(n_idx,:);
        memo_f = BSPTREE.f(n_idx,:);
        d = h;
        if memo_x(:, d) == x(idx, d) && flag == 0
            d1 = sum(xor(memo_x,BSPTREE.optimal));
            d2 = sum(xor(x(idx,:),BSPTREE.optimal));
            if d1 <= d2
                 memo_x(:, d) = ~memo_x(:, d);
            else
                 x(idx, d) =~x(idx, d);
            end
        else
            x(idx, d) =~x(idx, d);
            %--------------UNCODING------------------%
            y = zeros(D,1);
            m = inf*ones(D,Length);
            for i = 1:dim
                if mod(i,Length) == 0
                    m(i/Length,:) = memo_x(1,i-Length+1:i);
                end
            end
            for i = 1:D
                for j = 1:1:Length
                    y(i) = y(i) + m(i,j) * 2^(j-1);
                end
                xx(i) = (umax-umin)*y(i)/4095+umin;
            end
            %-----------END UNCODING-----------------%
            memo_f =  Test_Function(xx, D, testnumber);
        end
        BSPTREE.d(n_idx, 1) = d;
        
        l_idx = cur_bsptree_node + 1;
        r_idx = cur_bsptree_node + 2;
        cur_bsptree_node = cur_bsptree_node + 2;
        BSPTREE.l_idx(n_idx, 1) = l_idx;
        BSPTREE.r_idx(n_idx, 1) = r_idx;
        BSPTREE.p_idx(l_idx, 1) = n_idx;
        BSPTREE.p_idx(r_idx, 1) = n_idx;
        BSPTREE.d(l_idx, 1) = 0;
        BSPTREE.d(r_idx, 1) = 0;
        
        if x(idx, d) > 0.5
            BSPTREE.x(l_idx,: ) = x(idx,: );
            %--------------UNCODING------------------%
            y = zeros(D,1);
            m = inf*ones(D,Length);
            for i = 1:dim
                if mod(i,Length) == 0
                    m(i/Length,:) = x(idx,i-Length+1:i);
                end
            end
            for i = 1:D
                for j = 1:1:Length
                    y(i) = y(i) + m(i,j) * 2^(j-1);
                end
                xx(i) = (umax-umin)*y(i)/4095+umin;
            end
            %-----------END UNCODING-----------------%
            
            BSPTREE.f(l_idx,:) =  Test_Function(xx, D, testnumber);
            
            cur_popf(idx, 1) = BSPTREE.f(l_idx,: );
            if BSPTREE.f(l_idx,: )< BSPTREE.min_f
                BSPTREE.min_f = BSPTREE.f(l_idx,: );
                BSPTREE.optimial = BSPTREE.x(l_idx,: );
                BSPTREE.cro_d = d;
            end
            BSPTREE.x(r_idx,: ) = memo_x;
            BSPTREE.f(r_idx,: ) = memo_f;
        else
            BSPTREE.x(r_idx,: ) = x(idx,: );
            %--------------UNCODING------------------%
            y = zeros(D,1);
            m = inf*ones(D,Length);
            for i = 1:dim
                if mod(i,Length) == 0
                    m(i/Length,:) = x(idx,i-Length+1:i);
                end
            end
            for i = 1:D
                for j = 1:1:Length
                    y(i) = y(i) + m(i,j) * 2^(j-1);
                end
                xx(i) = (umax-umin)*y(i)/4095+umin;
            end
            %-----------END UNCODING-----------------%
            
            BSPTREE.f(r_idx,:) =  Test_Function(xx, D, testnumber);
            cur_popf(idx, 1) = BSPTREE.f(r_idx,: );
            if BSPTREE.f(r_idx,: )< BSPTREE.min_f
                BSPTREE.min_f = BSPTREE.f(r_idx,: );
                BSPTREE.optimial = BSPTREE.x(r_idx,: );
                BSPTREE.cro_d = d;
            end
            BSPTREE.x(l_idx,: ) = memo_x;
            BSPTREE.f(l_idx,: ) = memo_f;
        end
        flag = 0;
    end  % END OF IF-ELSE 'BSPTREE.D(1, 1) == -2' %

end  % END OF 'FOR IDX = 1 : POP_SIZE' %
cro_d = BSPTREE.cro_d;
end