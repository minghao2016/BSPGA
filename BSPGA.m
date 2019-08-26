clear;

global BSPTREE;
global test;
test = zeros(200,100);
% === USER-DEFINED PARAMETER DEFINITION ... === %

REPEAT_RUN      = 10; % REPEAT TIME %
MAX_FES         = 20000; % MAX EVALUATION TIME %
DIM             = 480;
save_tree_flag  = 'n';

% === NON-USER-DEFINED PARAMETER FOLLOWING ... === %

CR          = 0.5;  % THERE ARE ONLY TWOio PARAMETERS IN CNRGA % 
ProM1       = 0.05; % THE BSP TREE BASED LEARNING PROBABILITY %
ProM2       = 1/DIM; % MUTATION RATE %
POP_SIZE	= 100;  % THEY ARE: XOVER_RATE(CR) AND POP_SIZE  %
ROOT        = 1;    % THE ROOT NODE OF THE BSP TREE %

OUTPUT_FITNESS = Inf * ones(REPEAT_RUN, 1);  % FINAL RESULT OBTAINED %
CONVERGENCE_Gen = Inf * ones(REPEAT_RUN, 1); % COUNTER THE GEN OF CONVERGENCE %
TIME = Inf * ones(REPEAT_RUN, 1);
% === END OF PARAMETER DEFINITION === %

for r = 1 : REPEAT_RUN
     
    MAX_GENERATION   = floor(MAX_FES  / POP_SIZE);
    MAX_BSPTREE_NODE = 2 * MAX_FES;
      
% >> === === === INITIALIZATION === === === %

    cur_pop = round(rand(POP_SIZE, DIM));
    [BSPTREE, NUM_OF_BSPTREE] = bsptree_initilization(DIM, MAX_BSPTREE_NODE);
    [NUM_OF_BSPTREE, cur_pop, cur_fit, cross_d] = bsptree_insertion(cur_pop, NUM_OF_BSPTREE,1);
    T1 = cputime;
    bestff(1) = inf;
% >> === === === END OF INITIALIZATION === === === %

% >> === === === MAIN LOOP OF CNRGA === === === %
    for gen = 2 : MAX_GENERATION
        
        % --- --- --- UNIFORM CROSSOVER --- --- --- %
        i_idx = (1 : POP_SIZE)';
        rnd_idx = ceil(POP_SIZE * rand(POP_SIZE, 1));
        same_idx = find(rnd_idx(i_idx, 1) == i_idx);
        len_same_idx = length(same_idx);
        if (len_same_idx)
            for ls_idx = 1 : len_same_idx
                while rnd_idx(same_idx(ls_idx), 1) == i_idx(same_idx(ls_idx), 1)
                    rnd_idx(same_idx(ls_idx), 1) = ceil(POP_SIZE * rand(1, 1));
                end
            end
        end
        bin_eff = (rand(POP_SIZE, DIM) > CR);
        nxt_pop = cur_pop(i_idx,: ) .* bin_eff(:, :) + cur_pop(rnd_idx,: ) .* (~bin_eff(:, :));
        % --- --- --- END OF CROSSOVER --- --- --- %
        

        % --- --- --- MUTATION - 1--- --- --- --- --- %

        k = rand(POP_SIZE, 1);
        Temp = k<=ProM1;
        idex = find(Temp == 1);
        [~, cur_min_pop] = min(cur_fit);
        cross_region = cur_pop(cur_min_pop, 1:cross_d);
        [m1,n1] = size(cross_region);
        M1=size(idex,1);
        N1=1;
        midx = (1:m1)' * ones(1, M1);
        nidx = (1:n1)' * ones(1, N1);
        cross_regionb = cross_region(midx, nidx);
        nxt_pop(idex, 1:cross_d) = cross_regionb;

        % --- --- --- END OF MUTATION - 1 --- --- --- %
        
        % --- --- --- MUTATION - 2 --- --- --- --- ---%
%         PP = nxt_pop(:, 1:ceil(DIM/5));
%         k    = rand(POP_SIZE, ceil(DIM/5));
%         Temp = k<=ProM2;
%         PP(Temp) = 1-PP(Temp);
%         nxt_pop = [PP nxt_pop(:,ceil(DIM/5) + 1:DIM)];

        k    = rand(POP_SIZE, DIM);
        Temp = k<=ProM2;
        nxt_pop(Temp) = 1-nxt_pop(Temp);
        % --- --- --- END MUTATION - 2 --- --- --- ---%
        

        
        % --- --- --- INSERTION --- --- --- %
      
        [NUM_OF_BSPTREE, nxt_pop, nxt_fit, cross_d] = bsptree_insertion(nxt_pop, NUM_OF_BSPTREE,gen);
            
        
        % --- --- --- END: INSERTION--- --- %
        
        % --- --- --- SELECTION ----- %
        
        tmp_pop = [cur_pop; nxt_pop];
        tmp_fit = [cur_fit; nxt_fit];
        [~, sort_idx] = sort(tmp_fit);
        cur_pop(1 : POP_SIZE,: ) = tmp_pop(sort_idx(1 : POP_SIZE),: );
        cur_fit(:, 1) = tmp_fit(sort_idx(1 : POP_SIZE), 1);
        
        % --- --- --- END: SELECTION ----- %
        
        [cur_min_fit, cur_min_pop] = min(cur_fit);
        if cur_min_fit < OUTPUT_FITNESS(r, 1)
            OUTPUT_FITNESS(r, 1) = cur_min_fit;
            OUTPUT_POPULATION = cur_pop(cur_min_pop,:);
            CONVERGENCE_Gen(r, 1) = gen;
        end
        bestff(gen) = cur_min_fit;
    end  % END OF GEN = 2 : LIMIT_GENERATION %
    TIME(r, 1) = cputime - T1;
        plot(bestff);
    disp(['>> Fitness found (run ', int2str(r), '/', int2str(REPEAT_RUN), '): ', num2str(cur_min_fit)]);
    if strcmp(save_tree_flag, 'y') || strcmp(save_tree_flag, 'Y')
                fname = ['BSPTREE'];
        fname = strcat(fname);
        save(fname, 'BSPTREE');
        clear fname;
    end
end
MEANS  = mean(OUTPUT_FITNESS)
ERRORS = std(OUTPUT_FITNESS)
save_flag = 'n';
if strcmp(save_flag, 'y') || strcmp(save_flag, 'Y')
            fname = 'cNrGA';

    fname = strcat(fname);
    save(fname);
    clear save_flag fname;
end
    
