function [F_comm, F_sensing, feasible, SSNR_opt] = opt_jsc_GWO_variants(H_comm, sigmasq_comm, gamma, sensing_beamsteering, sensing_streams, sigmasq_sens, P_all, F_init, variant_name)
    SearchAgents_no = 30; 
    Max_iter = 50;       
    
    [U, M, N] = size(H_comm);
    num_vars = (U + sensing_streams) * (M * N); 
    dim = num_vars * 2; 
    
    Positions = rand(SearchAgents_no, dim) * 2 - 1; 
    
    if (strcmp(variant_name, 'HGWO') || strcmp(variant_name, 'IGWO') || strcmp(variant_name, 'WGWO')) && ~isempty(F_init)
        wolf_init_vector = encode_wolf(F_init, U, M, N, sensing_streams);
        Positions(1, :) = wolf_init_vector;
        for k = 2:5
             Positions(k, :) = wolf_init_vector + 0.1 * randn(1, dim);
        end
    end
    
    Alpha_pos = zeros(1, dim); Alpha_score = -inf; 
    Beta_pos = zeros(1, dim);  Beta_score = -inf;
    Delta_pos = zeros(1, dim); Delta_score = -inf;
    
    for l = 1:Max_iter
        if strcmp(variant_name, 'IGWO')
            a = 2 * (1 - (l / Max_iter)^2); 
        else
            a = 2 - l * (2 / Max_iter); 
        end
        
        for i = 1:SearchAgents_no
            Positions(i, :) = max(min(Positions(i, :), 10), -10);
            
            fitness = calculate_fitness(Positions(i, :), H_comm, sigmasq_comm, gamma, sensing_beamsteering, sigmasq_sens, P_all, U, M, N, sensing_streams);
            
            if fitness > Alpha_score
                Alpha_score = fitness; Alpha_pos = Positions(i, :);
            end
            if fitness > Beta_score && fitness < Alpha_score
                Beta_score = fitness; Beta_pos = Positions(i, :);
            end
            if fitness > Delta_score && fitness < Beta_score
                Delta_score = fitness; Delta_pos = Positions(i, :);
            end
        end
        
        for i = 1:SearchAgents_no
            for j = 1:dim
                r1 = rand(); r2 = rand();
                A1 = 2*a*r1 - a; C1 = 2*r2;
                D_alpha = abs(C1*Alpha_pos(j) - Positions(i,j));
                X1 = Alpha_pos(j) - A1*D_alpha;
                
                r1 = rand(); r2 = rand();
                A2 = 2*a*r1 - a; C2 = 2*r2;
                D_beta = abs(C2*Beta_pos(j) - Positions(i,j));
                X2 = Beta_pos(j) - A2*D_beta;
                
                r1 = rand(); r2 = rand();
                A3 = 2*a*r1 - a; C3 = 2*r2;
                D_delta = abs(C3*Delta_pos(j) - Positions(i,j));
                X3 = Delta_pos(j) - A3*D_delta;
                
                if strcmp(variant_name, 'WGWO')
                    Positions(i,j) = 0.5*X1 + 0.3*X2 + 0.2*X3;
                else
                    Positions(i,j) = (X1 + X2 + X3) / 3;
                end
            end
        end
    end
    
    [F_comm, F_sensing] = decode_wolf(Alpha_pos, U, M, N, sensing_streams);
    SSNR_opt = Alpha_score; 
    feasible = true;
end

function fit = calculate_fitness(pos, H_comm, sigmasq_comm, gamma, sensing_beamsteering, sigmasq_sens, P_all, U, M, N, streams)
    [F_comm, F_sensing] = decode_wolf(pos, U, M, N, streams);
    metrics = compute_metrics(H_comm, F_comm, sigmasq_comm, sensing_beamsteering, F_sensing, sigmasq_sens);
    
    obj = metrics.SSNR;
    penalty = 0;
    
    min_SINR = min(metrics.SINR, [], 'all'); 
    if min_SINR < gamma
        penalty = penalty + 1e6 * (gamma - min_SINR)^2;
    end
    
    total_power = sum(metrics.power, 'all');
    if total_power > P_all
        penalty = penalty + 1e6 * (total_power - P_all)^2;
    end
    
    fit = obj - penalty;
end

function [F_comm, F_sensing] = decode_wolf(pos, U, M, N, streams)
    len_vars = length(pos) / 2;
    pos_complex = pos(1:len_vars) + 1j * pos(len_vars+1:end);
    F_matrix_flat = reshape(pos_complex, M * N, U + streams);
    
    F_comm_flat = F_matrix_flat(:, 1:U);
    F_sensing_flat = F_matrix_flat(:, U+1:end);
    
    F_comm = permute(reshape(F_comm_flat, M, N, U), [3, 1, 2]);
    F_sensing = permute(reshape(F_sensing_flat, M, N, streams), [3, 1, 2]);
end

function vec = encode_wolf(F_matrix_3D, U, M, N, streams)
    F_flat_rows = reshape(F_matrix_3D, U+streams, M*N);
    F_flat_cols = F_flat_rows.'; 
    F_vec_complex = F_flat_cols(:);
    vec = [real(F_vec_complex); imag(F_vec_complex)].';
end
