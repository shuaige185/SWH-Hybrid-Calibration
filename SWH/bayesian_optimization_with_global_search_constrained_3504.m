% ===================== 针对3504组数据的约束增强多算法融合全局优化算法 =====================
function [best_params, para_range_10, history, optimization_info, all_method_params] = ...
    bayesian_optimization_with_global_search_constrained_3504(shuttle_calibration, Qs, max_iterations, constraint_weight)
    
    % 设置默认参数
    if nargin < 3
        max_iterations = 250;  
    end
    if nargin < 4
        constraint_weight = 0;  
    end
    
    % 参数边界
    nopt = 6;
    lb = [1.0000, 1.0000, 1.0000, 1.0000, 0.0010, 0.4500];
    ub = [5.0000, 5.0000, 1000.0000, 100.0000, 0.1000, 0.7000];
    
    best_params = [];  
    para_range_10 = [];  
    history = struct();  
    optimization_info = struct();  
    all_method_params = struct();  
    history.time = struct();  % 存储各算法时间
    
    % fmincon选项优化
    fmincon_options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point', ...
        'MaxIterations', 300,...         
        'MaxFunctionEvaluations', 5000, ...  
        'StepTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6, ...
        'FunctionTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-6, ...
        'UseParallel', false);        
    
    % 初始化历史记录和优化信息
    history.all_methods = struct();
    history.best_params_all = {};
    history.method_names = {};
    history.constraint_weight = constraint_weight;
    history.max_iterations = max_iterations;
    
    optimization_info.stages = {};
    optimization_info.best_method = '';
    optimization_info.all_results = struct();
    
    fprintf('===============================================\n');
    fprintf('针对3504组数据的约束增强多算法融合全局优化算法启动\n');
    fprintf('训练数据量: %d\n', size(shuttle_calibration, 1));
    fprintf('参数空间维度: %d\n', nopt);
    fprintf('最大迭代次数: %d\n', max_iterations);
    fprintf('===============================================\n\n');
    
    % ===================== 第一阶段：约束增强全局优化算法并行运行 =====================
    disp('第一阶段：运行多种约束增强的全局优化算法...');
    
    % 1. 约束增强贝叶斯优化
    disp('1. 约束增强贝叶斯优化...');
    t_start = tic; % ：记录时间
    [bayes_params, bayes_results] = run_constrained_bayesian_optimization_3504(...
        shuttle_calibration, Qs, 200, lb, ub, constraint_weight);  
    t_bayes = toc(t_start);
    fprintf('  => BO算法耗时: %.2f 秒\n', t_bayes);
     history.time.bayes = t_bayes;  
    history.all_methods.bayes = bayes_results;
    history.best_params_all{end+1} = bayes_params;
    history.method_names{end+1} = 'bayes';
    
    % 2. 约束增强遗传算法
    disp('2. 约束增强遗传算法...');
    t_start = tic;
    [ga_params, ga_results] = run_constrained_genetic_algorithm_optimization_3504(...
        shuttle_calibration, Qs, lb, ub, min(250, max_iterations), constraint_weight);
    t_ga = toc(t_start);
    fprintf('  => GA算法耗时: %.2f 秒\n', t_ga);
     history.time.ga = t_ga;  
    history.all_methods.ga = ga_results;
    history.best_params_all{end+1} = ga_params;
    history.method_names{end+1} = 'ga';
    
    % 3. 约束增强粒子群算法
    disp('3. 约束增强粒子群算法...');
    t_start = tic;
    [pso_params, pso_results] = run_constrained_particle_swarm_optimization_3504(...
        shuttle_calibration, Qs, lb, ub,min(250, max_iterations), constraint_weight);
    t_pso = toc(t_start);
    fprintf('  => PSO算法耗时: %.2f 秒\n', t_pso);
    history.time.pso = t_pso;
    history.all_methods.pso = pso_results;
    history.best_params_all{end+1} = pso_params;
    history.method_names{end+1} = 'pso';
    
    % 4. 约束增强麻雀算法
    disp('4. 约束增强麻雀算法...');
    t_start = tic;
    [ssa_params, ssa_results] = run_constrained_sparrow_search_algorithm_3504(...
        shuttle_calibration, Qs, lb, ub, min(250, max_iterations), constraint_weight);
    t_ssa = toc(t_start);
    fprintf('  => SSA算法耗时: %.2f 秒\n', t_ssa);
    history.time.ssa = t_ssa;
    history.all_methods.ssa = ssa_results;
    history.best_params_all{end+1} = ssa_params;
    history.method_names{end+1} = 'ssa';
    
    % ===================== 第二阶段：收集高质量参数作为fmincon起点 =====================
    disp('第二阶段：收集高质量参数作为fmincon起点...');
    
    all_top_params = {};
    all_top_fvals = {};
    
    for i = 1:length(history.method_names)
        method_name = history.method_names{i};
        method_results = history.all_methods.(method_name);
        
        if isfield(method_results, 'all_params') && isfield(method_results, 'all_fvals')
            % 去除异常高的失败结果
            good_indices = method_results.all_fvals < 1e7;
            good_params = method_results.all_params(good_indices, :);
            good_fvals = method_results.all_fvals(good_indices);
            
            % 用户要求：不足100个绝不随机补充，有多少用多少
            if size(good_params, 1) < 100
                fprintf('  %s: 警告 - 高质量参数不足100个，实际提取到 %d 个！\n', method_name, size(good_params, 1));
                top_params = good_params;
                top_fvals = good_fvals;
            else
                [sorted_fvals, sorted_idx] = sort(good_fvals);
                top_indices = sorted_idx(1:100);
                top_params = good_params(top_indices, :);
                top_fvals = sorted_fvals(1:100);
            end
            
            if isempty(top_params)
                error('致命错误: %s 方法无法提取任何高质量参数，程序终止！', method_name);
            end
            
            all_top_params{end+1} = top_params;
            all_top_fvals{end+1} = top_fvals;
            fprintf('  %s: 最终提取了 %d 个高质量参数作为起点\n', method_name, size(top_params, 1));
        else
            error('致命错误: %s 方法缺少历史参数输出，程序终止！', method_name);
        end
    end
    
    disp('对每类算法的参数运行 fmincon 局部微调...');
    all_fmincon_results = cell(length(history.method_names), 1);
    
    for i = 1:length(history.method_names)
        method_name = history.method_names{i};
        top_params = all_top_params{i};
        
        if ~isempty(top_params)
            disp(['  对 ' method_name ' 的参数运行fmincon...']);
            
            t_fmincon_start = tic; % 【新增】：开始记录 fmincon 耗时
            
            [fmincon_params, fmincon_fval, fmincon_info] = ...
                run_constrained_multi_start_fmincon_for_method_3504(...
                    top_params, shuttle_calibration, Qs, lb, ub, fmincon_options, constraint_weight);
                    
            t_fmincon_cost = toc(t_fmincon_start); % 【新增】：结束并计算 fmincon 耗时
            
            all_fmincon_results{i} = struct(...
                'method', method_name, ...
                'params', fmincon_params, ...
                'fval', fmincon_fval, ...
                'info', fmincon_info, ...
                'all_starting_points', top_params, ...        
                'all_params', fmincon_info.all_params, ...    
                'all_fvals', fmincon_info.all_fvals, ...
                'time_cost', t_fmincon_cost);        
            
            % 【修改】：在终端打印时，同时输出耗时
            fprintf('    %s-fmincon: 使用%d个起点微调完成. 最佳目标函数值 = %.4f | 微调耗时: %.2f 秒\n', ...
                method_name, size(top_params, 1), fmincon_fval, t_fmincon_cost);
            history.time.fmincon.(method_name) = t_fmincon_cost;
        end
    end
    
    % ===================== 第三阶段：融合结果，计算统计指标 =====================
  disp('第三阶段：融合原始结果和fmincon优化结果，计算Top10统计指标...');
    all_method_metrics = struct();
    
    for i = 1:length(history.method_names)
        method_name = history.method_names{i};
        method_results = history.all_methods.(method_name);
        fmincon_result = all_fmincon_results{i};
        
        original_params = method_results.all_params;
        original_fvals = method_results.all_fvals;
        
        fmincon_optimized_params = [];
        fmincon_optimized_fvals = [];
        
        if ~isempty(fmincon_result) && isfield(fmincon_result, 'all_params') && ~isempty(fmincon_result.all_params)
            fmincon_optimized_params = fmincon_result.all_params;
            fmincon_optimized_fvals = fmincon_result.all_fvals;
        end
        
        all_combined_params = [original_params; fmincon_optimized_params];
        all_combined_fvals = [original_fvals; fmincon_optimized_fvals];
        
        % 选择Top 10
        if ~isempty(all_combined_fvals)
            [sorted_fvals, sorted_idx] = sort(all_combined_fvals);
            n_top = min(10, length(sorted_fvals));
            top_idx = sorted_idx(1:n_top);
            top_params = all_combined_params(top_idx, :);
            top_fvals = sorted_fvals(1:n_top);
            
            is_fmincon_optimized = false(n_top, 1);
            if ~isempty(fmincon_optimized_params)
                for j = 1:n_top
                    param = top_params(j, :);
                    for k = 1:size(fmincon_optimized_params, 1)
                        if isequal(param, fmincon_optimized_params(k, :))
                            is_fmincon_optimized(j) = true;
                            break;
                        end
                    end
                end
            end
            
            [metrics_mean, metrics_all] = calculate_statistical_metrics_for_params_3504(...
                top_params, shuttle_calibration, Qs);
            
            all_method_metrics.(method_name) = struct(...
                'top10_params', top_params, ...
                'top10_fvals', top_fvals, ...
                'is_fmincon_optimized', is_fmincon_optimized, ...
                'metrics_mean', metrics_mean, ...
                'metrics_all', metrics_all, ...
                'num_fmincon_in_top10', sum(is_fmincon_optimized));
            
            fprintf('\n%s 方法 Top 10 参数统计 (包含%d个fmincon优化参数):\n', method_name, sum(is_fmincon_optimized));
            fprintf('  平均斜率: ET=%.4f\n', metrics_mean.slope_ET);
            fprintf('  平均R2: ET=%.4f\n', metrics_mean.R2_ET);
            fprintf('  平均RMSE: ET=%.4f\n', metrics_mean.RMSE_ET);
            fprintf('  平均NSE: ET=%.4f\n', metrics_mean.NSE_ET);
            fprintf('  平均MAE: ET=%.4f\n', metrics_mean.MAE_ET); 
        end
    end
    
    % ===================== 第四阶段：综合评估与选择最佳参数 =====================
    disp('第四阶段：综合评估与选择最佳参数...');
    all_candidates = [];
    
    for i = 1:length(all_fmincon_results)
        method_name = all_fmincon_results{i}.method;
        method_result = all_fmincon_results{i};
        if ~isempty(method_result.params)
            fval = unified_constrained_objective_3504(method_result.params, shuttle_calibration, Qs, constraint_weight);
            all_candidates = [all_candidates; struct('method', method_name, 'params', method_result.params, 'fval', fval, 'type', 'fmincon_optimized')];
        end
    end
    
    method_names = fieldnames(all_method_metrics);
    for i = 1:length(method_names)
        method_name = method_names{i};
        method_metrics = all_method_metrics.(method_name);
        if isfield(method_metrics, 'top10_params') && ~isempty(method_metrics.top10_params)
            mean_params = mean(method_metrics.top10_params, 1);
            fval = unified_constrained_objective_3504(mean_params, shuttle_calibration, Qs, constraint_weight);
            all_candidates = [all_candidates; struct('method', method_name, 'params', mean_params, 'fval', fval, 'type', 'mean_top10')];
        end
    end
    
    % 用户要求：仅按目标函数值进行纯净排序（越小越好）
    fvals = [all_candidates.fval];
    [~, sorted_idx] = sort(fvals);
    best_idx = sorted_idx(1);
    
    final_params = all_candidates(best_idx).params;
    best_fval = all_candidates(best_idx).fval;
    best_method = all_candidates(best_idx).method;
    best_type = all_candidates(best_idx).type;
    
    optimization_info.best_method = sprintf('%s (%s)', best_method, best_type);
    fprintf('\n? 选择 %s 方法的参数作为全场最佳参数 (目标函数值: %.6f)\n', optimization_info.best_method, best_fval);
    
    history.best_x = final_params;
    history.best_f = best_fval;
    history.all_candidates = all_candidates;
    history.all_method_metrics = all_method_metrics;
    
    [is_valid_final, final_validation] = validate_parameters_with_all_metrics_3504(final_params, shuttle_calibration, Qs);
    history.final_validation = final_validation;
    history.is_valid = is_valid_final;
    
    % ===================== 第五阶段：生成10组最优参数用于不确定性分析 =====================
   disp('第五阶段：提取最优参数用于不确定性分析...');
    all_best_params = [];
    for i = 1:length(method_names)
        method_name = method_names{i};
        method_metrics = all_method_metrics.(method_name);
        if isfield(method_metrics, 'top10_params') && ~isempty(method_metrics.top10_params)
            all_best_params = [all_best_params; method_metrics.top10_params];
        end
    end
    
    all_fvals = zeros(size(all_best_params, 1), 1);
    for i = 1:size(all_best_params, 1)
        all_fvals(i) = unified_constrained_objective_3504(all_best_params(i,:), shuttle_calibration, Qs, constraint_weight);
    end
    [~, sorted_idx] = sort(all_fvals);
    
    % 用户要求：不足10组绝不随机补充
    n_uncertainty = min(10, length(sorted_idx));
    if n_uncertainty < 10
        fprintf('  警告：收集的最优参数不足10组，实际只有 %d 组，将直接使用这些真实数据进行分析。\n', n_uncertainty);
    end
    
    top10_idx = sorted_idx(1:n_uncertainty);
    para_range_10_params = all_best_params(top10_idx, :);
    
    para_min = min(para_range_10_params, [], 1);
    para_max = max(para_range_10_params, [], 1);
    para_mean = mean(para_range_10_params, 1);
    para_sd = std(para_range_10_params, 0, 1);
    para_range_10 = [para_min; para_max; para_mean; para_sd];
    
    % ===================== 第六阶段：输出最终结果 =====================
   fprintf('\n===============================================\n');
    fprintf('针对3504组数据的多算法融合全局优化完成！\n');
    fprintf('训练数据量: %d\n', size(shuttle_calibration, 1));
    fprintf('最佳方法: %s\n', optimization_info.best_method);
    fprintf('最终目标函数值 (1-KGE距离, 越小越好): %.6f\n', best_fval);
    fprintf('最终参数: \n');
    fprintf('  b1: %.4f, b2: %.4f, b3: %.4f\n', final_params(1), final_params(2), final_params(3));
    fprintf('  a1: %.4f, g01: %.6f, LightExtCoef: %.4f\n', final_params(4), final_params(5), final_params(6));
    
   fprintf('\n验证结果:\n');
    fprintf('  ET斜率: %.4f\n', final_validation.slope_ET);
    fprintf('  ET-R2: %.4f\n', final_validation.R2_ET);
    fprintf('  ET-RMSE: %.4f\n', final_validation.RMSE_ET);
    fprintf('  ET-NSE: %.4f\n', final_validation.NSE_ET);
    fprintf('  ET-MAE: %.4f\n', final_validation.MAE_ET); 
    % ===================== 第七阶段：提取每个算法的参数集并对比【回应审稿人】 =====================
   disp('第七阶段：提取参数集并执行策略对比 (参数平均 vs 输出集合平均)...');
    method_list = {'bayes', 'ga', 'pso', 'ssa'};
    
    for m_idx = 1:length(method_list)
        method = method_list{m_idx};
        original_params = history.all_methods.(method).all_params;
        original_fvals = history.all_methods.(method).all_fvals;
        
        [~, sorted_idx] = sort(original_fvals);
        original_top100 = original_params(sorted_idx(1:min(100, length(sorted_idx))), :);
        
        fmincon_params = [];
        for i = 1:length(all_fmincon_results)
            if strcmp(all_fmincon_results{i}.method, method)
                fmincon_params = all_fmincon_results{i}.all_params; break;
            end
        end
        
        % 1. 原始Top10均值
        if size(original_top100, 1) >= 10
            top10_original = original_top100(1:10, :);
            param_set1 = mean(top10_original, 1);
            all_method_params.(method) = param_set1;
        end
        
        % 2. 混合Top10对比分析
        if ~isempty(fmincon_params)
            combined_params = [original_top100; fmincon_params];
            combined_fvals = zeros(size(combined_params, 1), 1);
            for i = 1:size(combined_params, 1)
                combined_fvals(i) = unified_constrained_objective_3504(combined_params(i,:), shuttle_calibration, Qs, constraint_weight);
            end
            [~, sorted_idx] = sort(combined_fvals);
            top10_combined = combined_params(sorted_idx(1:10), :);
            
            % 【策略B】参数算术平均 (Parameter Average)
            param_set2 = mean(top10_combined, 1);
            all_method_params.([method '_f']) = param_set2;
            [~, stats_avg] = validate_parameters_with_all_metrics_3504(param_set2, shuttle_calibration, Qs);
        fprintf('  %s_f (参数均值): R2=%.4f, NSE=%.4f, RMSE=%.4f, MAE=%.4f\n', method, stats_avg.R2_ET, stats_avg.NSE_ET, stats_avg.RMSE_ET, stats_avg.MAE_ET);
            
            % 【策略C】输出集合平均 (Output Ensemble) - 回应审稿人意见
           ET_sim_matrix = [];
            ET_obs_valid = [];
            for k = 1:10
                [shuttle_temp] = SWH_halfhour(shuttle_calibration, top10_combined(k,1), top10_combined(k,2), top10_combined(k,3), top10_combined(k,4), top10_combined(k,5), top10_combined(k,6), Qs);
                ET_sim = shuttle_temp(:,23);
                ET_obs = shuttle_temp(:,17);
                valid_idx = ET_obs ~= -99999;
                if isempty(ET_obs_valid), ET_obs_valid = ET_obs(valid_idx); end
                ET_sim_matrix = [ET_sim_matrix, ET_sim(valid_idx)];
            end
            ET_ens = mean(ET_sim_matrix, 2);
            stats_ens.RMSE_ET = sqrt(mean((ET_ens - ET_obs_valid).^2));
            stats_ens.MAE_ET = mean(abs(ET_ens - ET_obs_valid));
            stats_ens.NSE_ET = calculate_NSE(ET_ens, ET_obs_valid);
            coeCorr = corrcoef(ET_ens, ET_obs_valid);
            stats_ens.R2_ET = coeCorr(2,1)^2;
            fprintf('  %s_ens(输出集合平均): R2=%.4f, NSE=%.4f, RMSE=%.4f, MAE=%.4f\n\n', method, stats_ens.R2_ET, stats_ens.NSE_ET, stats_ens.RMSE_ET, stats_ens.MAE_ET);
        end
    end
    history.all_method_params = all_method_params;
    fprintf('======================================\n');
end

% ===================== 针对3504组数据的核心功能函数 =====================

function obj = unified_constrained_objective_3504(params, shuttle_data, Qs, ~)
    obj = objective_function_with_all_metrics_penalty_3504(params, shuttle_data, Qs);
end

% 真实目标函数 (基于 SGPI)
function obj = objective_function_with_all_metrics_penalty_3504(params, shuttle_data, Qs)
    lb = [1.0000, 1.0000, 1.0000, 1.0000, 0.0010, 0.4500];
    ub = [5.0000, 5.0000, 1000.0000, 100.0000, 0.1000, 0.7000];
    margin = 1e-4;
    
    % 严格的边界惩罚，防止参数越界
    if any(params < lb + margin) || any(params > ub - margin)
        obj = 100; 
        return;
    end
    
    try
        % 运行 SWH 模型
        [shuttle] = SWH_halfhour(shuttle_data, params(1), params(2), params(3), params(4), params(5), params(6), Qs);
        
        ET_sim = shuttle(:,23); 
        ET_obs = shuttle(:,17);
        
        % 数据清洗与匹配
        valid_idx = ET_obs ~= -99999 & ~isnan(ET_sim) & ~isinf(ET_sim);
        ET_sim = ET_sim(valid_idx); 
        ET_obs = ET_obs(valid_idx);
        
        % 若有效数据过少，施加严重惩罚
        if length(ET_sim) < 50
            obj = 100; 
            return; 
        end
        
        % 1. 计算 R2 (注意处理负相关或异常情况)
        r_mat = corrcoef(ET_sim, ET_obs); 
        if size(r_mat, 1) > 1 && ~isnan(r_mat(2,1))
            r2 = r_mat(2,1)^2; 
        else
            r2 = 0; 
        end
        
        % 2. 计算 NSE
        mean_obs_nse = mean(ET_obs);
        numerator = sum((ET_sim - ET_obs).^2);
        denominator = sum((ET_obs - mean_obs_nse).^2);
        if denominator == 0
            nse = -999;
        else
            nse = 1 - (numerator / denominator);
        end
        
        % 3. 计算相对 RMSE 和 MAE
        mean_obs = mean(ET_obs); 
        if mean_obs == 0, mean_obs = eps; end % 防止除以 0
        
        rmse = sqrt(mean((ET_sim - ET_obs).^2));
        mae = mean(abs(ET_sim - ET_obs));
        
        % 4. 构建 SGPI 代理目标函数 (越小越好)
        % 包含四个维度：(1-R2) + (1-NSE) + 相对RMSE + 相对MAE
        obj = (1 - r2) + (1 - nse) + (rmse / mean_obs) + (mae / mean_obs);
        
        if isnan(obj) || isinf(obj) || obj > 100
            obj = 100;
        end
    catch
        obj = 100;
    end
end
% 带所有指标验证的参数验证函数（针对3504组数据）
function [is_valid, results] = validate_parameters_with_all_metrics_3504(params, shuttle_data, Qs)
    try
        [shuttle_temp] = SWH_halfhour(shuttle_data, params(1), params(2), params(3), params(4), params(5), params(6), Qs);
        ET_sim = shuttle_temp(:,23);
        ET_obs = shuttle_temp(:,17);
        valid_idx = ET_obs ~= -99999;
        ET_sim = ET_sim(valid_idx);
        ET_obs = ET_obs(valid_idx);
        
        p_ET = polyfit(ET_obs, ET_sim, 1);  
        coeCorr_ET = corrcoef(ET_sim, ET_obs);
        if size(coeCorr_ET, 1) > 1, R_ET = coeCorr_ET(2,1)^2; else, R_ET = 0; end
        
        results.slope_ET = p_ET(1);
        results.R2_ET = R_ET;
        results.RMSE_ET = sqrt(mean((ET_sim - ET_obs).^2));
        results.MAE_ET = mean(abs(ET_sim - ET_obs)); % 新增 MAE
        results.NSE_ET = calculate_NSE(ET_sim, ET_obs);
        results.ET_samples = length(ET_sim);
        is_valid = true;
    catch
        is_valid = false;
        results.slope_ET = NaN; results.R2_ET = NaN; results.RMSE_ET = NaN; 
        results.NSE_ET = NaN; results.MAE_ET = NaN; results.ET_samples = 0;
    end
end

% 计算多组参数的统计指标矩阵
function [metrics_mean, metrics_all] = calculate_statistical_metrics_for_params_3504(params_matrix, shuttle_data, Qs)
    n_params = size(params_matrix, 1);
    metrics_all.slope_ET = zeros(n_params, 1);
    metrics_all.R2_ET = zeros(n_params, 1);
    metrics_all.RMSE_ET = zeros(n_params, 1);
    metrics_all.NSE_ET = zeros(n_params, 1);
    metrics_all.MAE_ET = zeros(n_params, 1); 
    
    for i = 1:n_params
        [~, stats] = validate_parameters_with_all_metrics_3504(params_matrix(i,:), shuttle_data, Qs);
        metrics_all.slope_ET(i) = stats.slope_ET;
        metrics_all.R2_ET(i) = stats.R2_ET;
        metrics_all.RMSE_ET(i) = stats.RMSE_ET;
        metrics_all.NSE_ET(i) = stats.NSE_ET;
        metrics_all.MAE_ET(i) = stats.MAE_ET;
    end
    
    metrics_mean.slope_ET = nanmean(metrics_all.slope_ET);
    metrics_mean.R2_ET = nanmean(metrics_all.R2_ET);
    metrics_mean.RMSE_ET = nanmean(metrics_all.RMSE_ET);
    metrics_mean.NSE_ET = nanmean(metrics_all.NSE_ET);
    metrics_mean.MAE_ET = nanmean(metrics_all.MAE_ET);
end

% 针对3504组数据的约束增强贝叶斯优化
function [best_params, results] = run_constrained_bayesian_optimization_3504(shuttle_data, Qs, max_evaluations, lb, ub, constraint_weight)
    fprintf('针对3504组数据的约束增强贝叶斯优化开始...\n');
    
    try
        % 准备贝叶斯优化变量
        optimVars = [
            optimizableVariable('b1', [lb(1), ub(1)], 'Type', 'real'),
            optimizableVariable('b2', [lb(2), ub(2)], 'Type', 'real'),
            optimizableVariable('b3', [lb(3), ub(3)], 'Type', 'real'),
            optimizableVariable('a1', [lb(4), ub(4)], 'Type', 'real'),
            optimizableVariable('g01', [lb(5), ub(5)], 'Type', 'real'),
            optimizableVariable('LightExtCoef', [lb(6), ub(6)], 'Type', 'real')
        ];
        
        % 定义约束增强的贝叶斯优化目标函数
        bayes_constrained_objective = @(params) constrained_bayes_objective_wrapper_3504(params, shuttle_data, Qs, constraint_weight);
        
        % 运行贝叶斯优化（针对3504组数据增加评估次数）
        results_bayes = bayesopt(bayes_constrained_objective, optimVars, ...
            'MaxObjectiveEvaluations', max_evaluations, ...
            'AcquisitionFunctionName', 'expected-improvement-plus', ...
            'IsObjectiveDeterministic', true, ...
            'ExplorationRatio', 0.6, ...
            'NumSeedPoints', 15, ...
            'Verbose', 1, ...
            'PlotFcn', []);
        
        % 获取结果
        best_params = table2array(results_bayes.XAtMinObjective);
        best_fval = results_bayes.MinObjective;
        
        % 提取所有历史点
        all_params = table2array(results_bayes.XTrace);
        all_fvals = results_bayes.ObjectiveTrace;
        
        % 验证最佳参数是否满足约束
        [c_best, ~] = parameter_constraints_with_all_metrics_3504(best_params, shuttle_data, Qs);
        if any(c_best > 0)
            fprintf('警告：贝叶斯优化最佳参数违反约束\n');
            
            % 从历史点中找到满足约束的最佳参数
            feasible_indices = [];
            for i = 1:size(all_params, 1)
                [c_temp, ~] = parameter_constraints_with_all_metrics_3504(all_params(i,:), shuttle_data, Qs);
                if ~any(c_temp > 0)
                    feasible_indices = [feasible_indices; i];
                end
            end
            
            if ~isempty(feasible_indices)
                feasible_fvals = all_fvals(feasible_indices);
                [best_feasible_fval, best_idx] = min(feasible_fvals);
                best_params = all_params(feasible_indices(best_idx), :);
                best_fval = best_feasible_fval;
                fprintf('已选择满足约束的次优解，目标函数值: %.4f\n', best_fval);
            end
        end
        
    catch ME
        fprintf('约束增强贝叶斯优化失败: %s，使用简化版本替代\n', ME.message);
        
        % 使用带约束的随机搜索替代
        n_samples = max_evaluations;
        all_params = zeros(n_samples, length(lb));
        all_fvals = zeros(n_samples, 1);
        
        best_fval = Inf;
        best_params = [];
        
        for i = 1:n_samples
            params = zeros(1, length(lb));
            for j = 1:length(lb)
                params(j) = lb(j) + (ub(j) - lb(j)) * rand();
            end
            
            fval = unified_constrained_objective_3504(params, shuttle_data, Qs, constraint_weight);
            
            all_params(i, :) = params;
            all_fvals(i) = fval;
            
            if fval < best_fval
                best_fval = fval;
                best_params = params;
            end
        end
    end
    
    % 保存结果
    results = struct();
    results.best_fval = best_fval;
    results.all_params = all_params;
    results.all_fvals = all_fvals;
    results.iterations = max_evaluations;
    
    fprintf('约束增强贝叶斯优化完成！最佳目标函数值: %f\n', best_fval);
end

% 贝叶斯优化包装器（针对3504组数据）
function obj = constrained_bayes_objective_wrapper_3504(params_table, shuttle_data, Qs, constraint_weight)
    try
        % 将表格转换为向量
        if istable(params_table)
            params = [params_table.b1, params_table.b2, params_table.b3, ...
                      params_table.a1, params_table.g01, params_table.LightExtCoef];
        else
            params = params_table;
        end
        
        % 调用统一约束增强目标函数
        obj = unified_constrained_objective_3504(params, shuttle_data, Qs, constraint_weight);
        
        if ~isscalar(obj)
            obj = 1e10;
        end
        
    catch
        obj = 1e10;
    end
end

% 针对3504组数据的约束增强遗传算法（增加种群大小）
function [best_params, results] = run_constrained_genetic_algorithm_optimization_3504(shuttle_data, Qs, lb, ub, max_iterations, constraint_weight)
    fprintf('针对3504组数据的约束增强遗传算法优化开始...\n');
    
    % ========== 增强的算法参数 ==========
    n_pop = 200;           % 种群大小
    n_dim = length(lb);
    n_elite = 20;          % 增加精英数量
    eta_c = 20;            % 交叉分布指数（用于SBX）
    eta_m = 20;            % 变异分布指数（用于多项式变异）
    mutation_rate_initial = 0.3;  % 初始变异率
    crossover_rate_initial = 0.9; % 初始交叉率
    
    % ========== 添加：自适应参数公式 ==========
    % 交叉率和变异率随迭代自适应变化
    crossover_rate = @(iter) crossover_rate_initial - ...
        (crossover_rate_initial - 0.6) * (iter / max_iterations);
    
    mutation_rate = @(iter) mutation_rate_initial * ...
        exp(-5 * iter / max_iterations);
    
    % ========== 初始化种群 ==========
    population = zeros(n_pop, n_dim);
    fitness = zeros(n_pop, 1);
    feasibility = false(n_pop, 1);
    
    for i = 1:n_pop
        for j = 1:n_dim
            population(i, j) = lb(j) + (ub(j) - lb(j)) * rand();
        end
        fitness(i) = unified_constrained_objective_3504(population(i,:), shuttle_data, Qs, constraint_weight);
        [c_temp, ~] = parameter_constraints_with_all_metrics_3504(population(i,:), shuttle_data, Qs);
        feasibility(i) = ~any(c_temp > 0);
    end
    
    all_params_history = population;
    all_fvals_history = fitness;
    feasibility_history = zeros(max_iterations, 1);
    
    % 初始化停滞检测变量
    best_fitness_prev = min_fitness_with_feasibility(fitness, feasibility);
    stall_counter = 0;
    stall_threshold = 1e-4;
    stall_limit = 10;
    
    % ========== 迭代优化 ==========
    for iter = 1:max_iterations
        % 获取当前代的自适应参数
        current_crossover_rate = crossover_rate(iter);
        current_mutation_rate = mutation_rate(iter);
        
        % 统计可行性
        feasible_count = sum(feasibility);
        feasibility_history(iter) = feasible_count / n_pop;
        
        % ========== 增强的选择策略：排名选择 ==========
        [sorted_fitness, sorted_idx] = sort_fitness_with_feasibility(fitness, feasibility);
        rank = 1:n_pop;
        
        % 线性排名选择概率（选择压力s=1.5）
        s = 1.5;  % 选择压力参数
        max_rank = n_pop;
        min_rank = 1;
        selection_prob = (2 - s)/max_rank + 2*(s - 1)*(max_rank - rank + 1) / (max_rank*(max_rank - 1));
        selection_prob = selection_prob / sum(selection_prob);
        
        % 精英选择
        new_population = zeros(n_pop, n_dim);
        new_fitness = zeros(n_pop, 1);
        new_feasibility = false(n_pop, 1);
        
        % 保留精英（前n_elite个）
        for i = 1:n_elite
            new_population(i,:) = population(sorted_idx(i),:);
            new_fitness(i) = fitness(sorted_idx(i));
            new_feasibility(i) = feasibility(sorted_idx(i));
        end
        
        % ========== 交叉和变异产生新个体 ==========
        for i = (n_elite+1):2:n_pop
            % 选择父代（排名选择）
            parent1_idx = randsample(1:n_pop, 1, true, selection_prob);
            parent2_idx = randsample(1:n_pop, 1, true, selection_prob);
            
            parent1 = population(parent1_idx, :);
            parent2 = population(parent2_idx, :);
            
            % 初始化子代
            child1 = parent1;
            child2 = parent2;
            
            % ========== 增强的交叉：模拟二进制交叉（SBX） ==========
            if rand() < current_crossover_rate
                for j = 1:n_dim
                    if rand() <= 0.5
                        % SBX交叉
                        u = rand();
                        if u <= 0.5
                            beta = (2*u)^(1/(eta_c+1));
                        else
                            beta = (1/(2*(1-u)))^(1/(eta_c+1));
                        end
                        
                        child1(j) = 0.5 * ((1+beta)*parent1(j) + (1-beta)*parent2(j));
                        child2(j) = 0.5 * ((1-beta)*parent1(j) + (1+beta)*parent2(j));
                        
                        % 边界检查
                        child1(j) = max(lb(j), min(ub(j), child1(j)));
                        child2(j) = max(lb(j), min(ub(j), child2(j)));
                    end
                end
            end
            
            % ========== 增强的变异：多项式变异 ==========
            if rand() < current_mutation_rate
                for child_idx = 1:2
                    if child_idx == 1
                        child = child1;
                    else
                        child = child2;
                    end
                    
                    for j = 1:n_dim
                        if rand() < 1/n_dim  % 每个维度有一定概率变异
                            delta = 0;
                            u = rand();
                            if u <= 0.5
                                delta = (2*u)^(1/(eta_m+1)) - 1;
                            else
                                delta = 1 - (2*(1-u))^(1/(eta_m+1));
                            end
                            
                            child(j) = child(j) + delta * (ub(j) - lb(j));
                            child(j) = max(lb(j), min(ub(j), child(j)));
                        end
                    end
                    
                    if child_idx == 1
                        child1 = child;
                    else
                        child2 = child;
                    end
                end
            end
            
            % 评估子代
            new_population(i,:) = child1;
            new_fitness(i) = unified_constrained_objective_3504(child1, shuttle_data, Qs, constraint_weight);
            
            if i+1 <= n_pop
                new_population(i+1,:) = child2;
                new_fitness(i+1) = unified_constrained_objective_3504(child2, shuttle_data, Qs, constraint_weight);
            end
            
            % 检查可行性
            [c1, ~] = parameter_constraints_with_all_metrics_3504(child1, shuttle_data, Qs);
            [c2, ~] = parameter_constraints_with_all_metrics_3504(child2, shuttle_data, Qs);
            
            new_feasibility(i) = ~any(c1 > 0);
            if i+1 <= n_pop
                new_feasibility(i+1) = ~any(c2 > 0);
            end
        end
        
        % ========== 替换种群 ==========
        population = new_population;
        fitness = new_fitness;
        feasibility = new_feasibility;
        
        % 记录历史
        all_params_history = [all_params_history; population];
        all_fvals_history = [all_fvals_history; fitness];
        
        % 检测停滞
        current_best = min_fitness_with_feasibility(fitness, feasibility);
        if abs(current_best - best_fitness_prev) < stall_threshold
            stall_counter = stall_counter + 1;
        else
            stall_counter = 0;
            best_fitness_prev = current_best;
        end
        
        % 停滞时执行自适应调整
        if stall_counter >= stall_limit
            fprintf('  迭代 %d: 检测到停滞，执行自适应调整（变异率提高、重启部分个体）\n', iter);
            
            % 1. 提高变异率
            current_mutation_rate = min(0.5, current_mutation_rate * 1.5);
            
            % 2. 重新初始化后20%的最差个体（跳过精英）
            n_reinit = ceil(0.2 * n_pop);
            [~, worst_idx] = sort(fitness, 'descend');
            for k = 1:n_reinit
                idx = worst_idx(k);
                if idx > n_elite  % 不重新初始化精英
                    for j = 1:n_dim
                        population(idx, j) = lb(j) + (ub(j) - lb(j)) * rand();
                    end
                    fitness(idx) = unified_constrained_objective_3504(population(idx,:), shuttle_data, Qs, constraint_weight);
                    [c_temp, ~] = parameter_constraints_with_all_metrics_3504(population(idx,:), shuttle_data, Qs);
                    feasibility(idx) = ~any(c_temp > 0);
                end
            end
            
            stall_counter = 0;  % 重置计数器
        end
        
        % 显示进度
        if mod(iter, 20) == 0 || iter <= 5
            best_feasible_idx = find_best_feasible(fitness, feasibility);
            fprintf('  迭代 %d/%d, 交叉率: %.2f, 变异率: %.2f, 可行解比例: %.1f%%, 最佳适应度: %.4f\n', ...
                iter, max_iterations, current_crossover_rate, current_mutation_rate, ...
                feasibility_history(iter)*100, fitness(best_feasible_idx));
        end
    end
    
    % 找到最佳参数
    best_feasible_idx = find_best_feasible(fitness, feasibility);
    best_params = population(best_feasible_idx, :);
    best_fval = fitness(best_feasible_idx);
    
    % 保存结果
    results = struct();
    results.best_fval = best_fval;
    results.all_params = all_params_history;
    results.all_fvals = all_fvals_history;
    results.iterations = max_iterations;
    results.feasibility_history = feasibility_history;
    results.is_feasible = feasibility(best_feasible_idx);
    
    fprintf('增强的约束增强遗传算法优化完成！最佳目标函数值: %f, 可行: %s\n', ...
        best_fval, string(results.is_feasible));
end

% 针对3504组数据的约束增强粒子群算法（增加粒子数量）
function [best_params, results] = run_constrained_particle_swarm_optimization_3504(shuttle_data, Qs, lb, ub, max_iterations, constraint_weight)
    fprintf('针对3504组数据的增强约束增强粒子群优化算法开始...\n');
    
    % ========== 增强的算法参数 ==========
    n_particles = 200;
    n_dim = length(lb);
    
    % ========== 添加：动态邻域拓扑 ==========
    topology_type = 'ring';  % 'ring', 'star', 'fully_connected', 'dynamic'
    neighborhood_size = 5;   % 环形拓扑的邻居数量
    
    % ========== 初始化粒子 ==========
    particles = zeros(n_particles, n_dim);
    velocity = zeros(n_particles, n_dim);
    pbest = zeros(n_particles, n_dim);
    pbest_fval = zeros(n_particles, 1);
    pbest_feasible = false(n_particles, 1);
    
    % 初始化记录
    all_params_history = [];
    all_fvals_history = [];
    feasibility_history = zeros(max_iterations, 1);
    
    for i = 1:n_particles
        for j = 1:n_dim
            particles(i, j) = lb(j) + (ub(j) - lb(j)) * rand();
            velocity(i, j) = (ub(j) - lb(j)) * (rand() - 0.5) * 0.1;
        end
        
        pbest(i, :) = particles(i, :);
        pbest_fval(i) = unified_constrained_objective_3504(particles(i,:), shuttle_data, Qs, constraint_weight);
        [c_temp, ~] = parameter_constraints_with_all_metrics_3504(particles(i,:), shuttle_data, Qs);
        pbest_feasible(i) = ~any(c_temp > 0);
    end
    
    % ========== 找到全局最优 ==========
    [gbest_fval, gbest_idx] = min(pbest_fval);
    gbest = pbest(gbest_idx, :);
    gbest_feasible = pbest_feasible(gbest_idx);
    
    % 调整：如果没有可行解，选择目标函数值最小的
    if ~gbest_feasible
        feasible_indices = find(pbest_feasible);
        if ~isempty(feasible_indices)
            [gbest_fval, best_feasible_idx] = min(pbest_fval(feasible_indices));
            gbest = pbest(feasible_indices(best_feasible_idx), :);
            gbest_feasible = true;
            gbest_idx = feasible_indices(best_feasible_idx);
        end
    end
    
    % 记录历史
    all_params_history = [all_params_history; particles];
    all_fvals_history = [all_fvals_history; pbest_fval];
    
    % ========== 迭代优化 ==========
    for iter = 1:max_iterations
        % ========== 动态调整算法参数 ==========
        % 惯性权重（线性递减）
        w_max = 0.9;
        w_min = 0.4;
        w = w_max - (w_max - w_min) * (iter / max_iterations);
        
        % 学习因子（时变）
        c1_max = 2.5; c1_min = 1.5;
        c2_max = 2.5; c2_min = 1.5;
        c1 = c1_max - (c1_max - c1_min) * (iter / max_iterations);
        c2 = c2_min + (c2_max - c2_min) * (iter / max_iterations);
        
        % ========== 添加：收缩因子（Clerc & Kennedy, 2002） ==========
        phi = c1 + c2;
        if phi > 4
            kappa = 1;
            chi = 2 * kappa / abs(2 - phi - sqrt(phi^2 - 4*phi));
        else
            chi = 1;
        end
        
        % 统计可行性
        feasible_count = 0;
        
       % ========== 更新每个粒子 ==========（Kennedy, J., & Mendes, R. (2002). Population structure and particle swarm performance.Proceedings of the 2002 Congress on Evolutionary Computation.）
        for i = 1:n_particles
            % 获取邻居（动态拓扑）
            switch topology_type
                case 'ring'
                    % 环形拓扑
                    neighbors = mod((i-neighborhood_size):(i+neighborhood_size), n_particles);
                    neighbors(neighbors == 0) = n_particles;
                    neighbors_gbest = get_best_neighbor(pbest_fval, pbest_feasible, neighbors);
                    
                case 'star'
                    % 星形拓扑（所有粒子连接到gbest）
                    neighbors_gbest = gbest;
                    
                case 'fully_connected'
                    % 全连接（标准PSO）
                    neighbors_gbest = gbest;
                    
                case 'dynamic'
                    % 动态拓扑（每50代改变一次）
                    if mod(iter, 50) == 0
                        topology_type_list = {'ring', 'star', 'fully_connected'};
                        topology_type = topology_type_list{randi([1,3])};
                    end
                    neighbors_gbest = gbest;
            end
            
            % ========== 标准PSO速度更新（带收缩因子） ==========
            r1 = rand(1, n_dim);
            r2 = rand(1, n_dim);
            
            velocity(i,:) = chi * (w * velocity(i,:) + ...
                c1 * r1 .* (pbest(i,:) - particles(i,:)) + ...
                c2 * r2 .* (neighbors_gbest - particles(i,:)));
            
            % ========== 添加：速度限制 ==========
            v_max_factor = 0.15;
            v_max = (ub - lb) * v_max_factor;
            for j = 1:n_dim
                if velocity(i, j) > v_max(j)
                    velocity(i, j) = v_max(j);
                elseif velocity(i, j) < -v_max(j)
                    velocity(i, j) = -v_max(j);
                end
            end
            
            % 更新位置
            particles(i,:) = particles(i,:) + velocity(i,:);
            
            % 边界检查
            for j = 1:n_dim
                if particles(i, j) < lb(j)
                    particles(i, j) = lb(j);
                    velocity(i, j) = -0.5 * velocity(i, j);  % 边界反弹
                elseif particles(i, j) > ub(j)
                    particles(i, j) = ub(j);
                    velocity(i, j) = -0.5 * velocity(i, j);  % 边界反弹
                end
            end
            
            % 计算新适应度
            new_fval = unified_constrained_objective_3504(particles(i,:), shuttle_data, Qs, constraint_weight);
            
            % 检查可行性
            [c_new, ~] = parameter_constraints_with_all_metrics_3504(particles(i,:), shuttle_data, Qs);
            is_feasible = ~any(c_new > 0);
            if is_feasible
                feasible_count = feasible_count + 1;
            end
            
            % ========== 更新个体最优 ==========
            update_pbest = false;
            if is_feasible && ~pbest_feasible(i)
                update_pbest = true;
            elseif is_feasible && pbest_feasible(i)
                if new_fval < pbest_fval(i)
                    update_pbest = true;
                end
            elseif ~is_feasible && ~pbest_feasible(i)
                % 如果都不可行，选择约束违反更少的
                old_violation = sum(max(0, parameter_constraints_with_all_metrics_3504(pbest(i,:), shuttle_data, Qs)));
                new_violation = sum(max(0, c_new));
                if new_violation < old_violation || new_fval < pbest_fval(i)
                    update_pbest = true;
                end
            end
            
            if update_pbest
                pbest(i,:) = particles(i,:);
                pbest_fval(i) = new_fval;
                pbest_feasible(i) = is_feasible;
                
                % 更新全局最优
                update_gbest = false;
                if is_feasible && ~gbest_feasible
                    update_gbest = true;
                elseif is_feasible && gbest_feasible
                    if new_fval < gbest_fval
                        update_gbest = true;
                    end
                elseif ~is_feasible && ~gbest_feasible
                    if new_fval < gbest_fval
                        update_gbest = true;
                    end
                end
                
                if update_gbest
                    gbest = particles(i,:);
                    gbest_fval = new_fval;
                    gbest_feasible = is_feasible;
                    gbest_idx = i;
                end
            end
        end
        
        % ========== 添加：多样性保持机制 ==========
        if mod(iter, 10) == 0 && feasible_count < 0.1 * n_particles
            fprintf('  低多样性，重新初始化部分粒子...\n');
            % 重新初始化最差的20%粒子
            [~, worst_idx] = sort(pbest_fval, 'descend');
            reinit_num = ceil(0.2 * n_particles);
            
            for k = 1:reinit_num
                idx = worst_idx(k);
                for j = 1:n_dim
                    particles(idx, j) = lb(j) + (ub(j) - lb(j)) * rand();
                    velocity(idx, j) = (ub(j) - lb(j)) * (rand() - 0.5) * 0.1;
                end
            end
        end
        
        % 记录可行性历史
        feasibility_history(iter) = feasible_count / n_particles;
        
        % 记录历史
        all_params_history = [all_params_history; particles];
        all_fvals_history = [all_fvals_history; pbest_fval];
        
        % 显示进度
        if mod(iter, 20) == 0 || iter <= 5
            fprintf('  迭代 %d/%d, w=%.2f, c1=%.2f, c2=%.2f, 可行解比例: %.1f%%, 最佳适应度: %.4f\n', ...
                iter, max_iterations, w, c1, c2, feasibility_history(iter)*100, gbest_fval);
        end
    end
    
    % 返回最佳参数
    best_params = gbest;
    
    % 保存结果
    results = struct();
    results.best_fval = gbest_fval;
    results.all_params = all_params_history;
    results.all_fvals = all_fvals_history;
    results.iterations = max_iterations;
    results.feasibility_history = feasibility_history;
    results.is_feasible = gbest_feasible;
    
    fprintf('增强的约束增强粒子群优化算法完成！最佳目标函数值: %f, 可行: %s\n', ...
        gbest_fval, string(gbest_feasible));
end

% ========== 辅助函数：获取最佳邻居 ==========
function best_neighbor = get_best_neighbor(fvals, feasibility, neighbor_indices)
    % 从邻居中找到最佳个体（可行解优先）
    neighbor_fvals = fvals(neighbor_indices);
    neighbor_feasibility = feasibility(neighbor_indices);
    
    % 先找可行解
    feasible_neighbors = neighbor_indices(neighbor_feasibility);
    if ~isempty(feasible_neighbors)
        [~, idx] = min(fvals(feasible_neighbors));
        best_neighbor = feasible_neighbors(idx);
    else
        % 没有可行解，选择适应度最好的
        [~, idx] = min(neighbor_fvals);
        best_neighbor = neighbor_indices(idx);
    end
end

% 针对3504组数据的约束增强麻雀算法（增加种群大小）
function [best_params, results] = run_constrained_sparrow_search_algorithm_3504(shuttle_data, Qs, lb, ub, max_iterations, constraint_weight)
    fprintf('针对3504组数据的增强约束增强麻雀搜索算法开始...\n');
    
    % ========== 标准SSA参数 ==========
    n_pop = 200;
    n_dim = length(lb);
    
    % 标准SSA参数
    pd = 0.7;        % 发现者比例
    sd = 0.2;        % 警戒者比例
    ST = 0.6;        % 安全阈值
    R2_min = 0.8;    % 预警值最小值
    R2_max = 1.0;    % 预警值最大值
    
    alpha = 1.0;     % 发现者更新参数
    Q_max = 1.0;     % 随机数最大值
    Q_min = 0.0;     % 随机数最小值
    
    % ========== 初始化种群 ==========
    population = zeros(n_pop, n_dim);
    fitness = zeros(n_pop, 1);
    is_feasible = false(n_pop, 1);
    
    % 使用拉丁超立方采样（更好的初始化）
    population = lhsdesign(n_pop, n_dim);
    for i = 1:n_pop
        for j = 1:n_dim
            population(i, j) = lb(j) + (ub(j) - lb(j)) * population(i, j);
        end
        fitness(i) = unified_constrained_objective_3504(population(i,:), shuttle_data, Qs, constraint_weight);
        [c_temp, ~] = parameter_constraints_with_all_metrics_3504(population(i,:), shuttle_data, Qs);
        is_feasible(i) = ~any(c_temp > 0);
    end
    
    % 记录历史
    all_params_history = population;
    all_fvals_history = fitness;
    feasibility_history = zeros(max_iterations, 1);
    
    % ========== 标准SSA迭代 ==========
    for iter = 1:max_iterations
        % 动态调整预警值
        R2 = R2_min + (R2_max - R2_min) * rand();
        
        % 计算发现者数量和警戒者数量
        n_discoverers = round(pd * n_pop);
        n_watchers = round(sd * n_pop);
        
        % ========== 步骤1：对麻雀按适应度排序 ==========
        [sorted_fitness, sorted_idx] = sort_fitness_with_feasibility(fitness, is_feasible);
        best_fitness = sorted_fitness(1);
        worst_fitness = sorted_fitness(end);
        best_idx = sorted_idx(1);
        worst_idx = sorted_idx(end);
        
        best_solution = population(best_idx, :);
        worst_solution = population(worst_idx, :);
        
        % 统计可行性
        feasible_count = sum(is_feasible);
        feasibility_history(iter) = feasible_count / n_pop;
        
        % ========== 步骤2：发现者更新（标准SSA公式） ==========
        for i = 1:n_discoverers
            if R2 < ST
                % 安全状态，广泛搜索
                % 公式：X_{i,j}^{t+1} = X_{i,j}^t · exp(-i/(α·T_max))
                for j = 1:n_dim
                    population(i, j) = population(i, j) * exp(-i / (alpha * max_iterations));
                end
            else
                % 危险状态，飞向安全区域
                % 公式：X_{i,j}^{t+1} = X_{i,j}^t + Q·L
                L = ones(1, n_dim);
                Q = Q_min + (Q_max - Q_min) * randn();
                for j = 1:n_dim
                    population(i, j) = population(i, j) + Q * L(j);
                end
            end
            
            % 边界检查
            for j = 1:n_dim
                population(i, j) = max(lb(j), min(ub(j), population(i, j)));
            end
        end
        
        % ========== 步骤3：跟随者更新（标准SSA公式） ==========
        for i = (n_discoverers+1):n_pop
            if i > n_pop/2
                % 较差跟随者：饥饿状态
                % 公式：X_{i,j}^{t+1} = Q · exp((X_{worst} - X_{i,j}^t)/i?)
                Q = Q_min + (Q_max - Q_min) * randn();
                for j = 1:n_dim
                    population(i, j) = Q * exp((worst_solution(j) - population(i, j)) / i^2);
                end
            else
                % 较好跟随者：向最优位置移动 (数学简化版，极大提升速度)
                % 参考文献: Xue, J., & Shen, B. (2020). Sparrow search algorithm.
                A = randi([0, 1], 1, n_dim) * 2 - 1; 
                A_plus = A' / (A * A'); 
                L = ones(1, n_dim);
                difference = abs(population(i, :) - best_solution);
                update_term = difference .* (A_plus * L)';
                population(i, :) = best_solution + update_term(1, :);
            end
            
            % 边界检查
            for j = 1:n_dim
                population(i, j) = max(lb(j), min(ub(j), population(i, j)));
            end
        end
        
        % ========== 步骤4：警戒者更新（标准SSA公式） ==========
        % 随机选择警戒者
        watcher_indices = randperm(n_pop, n_watchers);
        
        for idx = 1:length(watcher_indices)
            i = watcher_indices(idx);
            
            % 计算当前麻雀的适应度
            current_fval = fitness(i);
            mean_fitness = mean(fitness);
            
            if current_fval > mean_fitness
                % 公式：X_{i,j}^{t+1} = X_{best}^t + β·|X_{i,j}^t - X_{best}^t|
                beta = randn();
                for j = 1:n_dim
                    population(i, j) = best_solution(j) + beta * abs(population(i, j) - best_solution(j));
                end
            else
                % 公式：X_{i,j}^{t+1} = X_{i,j}^t + K·(|X_{i,j}^t - X_{worst}^t|/(f_i - f_w + ε))
                K = 2 * rand() - 1;  % [-1, 1]
                epsilon = 1e-8;
                
                for j = 1:n_dim
                    numerator = K * abs(population(i, j) - worst_solution(j));
                    denominator = current_fval - worst_fitness + epsilon;
                    population(i, j) = population(i, j) + numerator / denominator;
                end
            end
            
            % 边界检查
            for j = 1:n_dim
                population(i, j) = max(lb(j), min(ub(j), population(i, j)));
            end
        end
        
        % ========== 评估新种群 ==========
        for i = 1:n_pop
            % 计算新适应度
            new_fval = unified_constrained_objective_3504(population(i,:), shuttle_data, Qs, constraint_weight);
            
            % 检查可行性
            [c_new, ~] = parameter_constraints_with_all_metrics_3504(population(i,:), shuttle_data, Qs);
            new_feasible = ~any(c_new > 0);
            
            % 更新个体
            update_individual = false;
            
            if new_feasible && ~is_feasible(i)
                update_individual = true;
            elseif new_feasible && is_feasible(i)
                if new_fval < fitness(i)
                    update_individual = true;
                end
            elseif ~new_feasible && ~is_feasible(i)
                old_violation = sum(max(0, parameter_constraints_with_all_metrics_3504(population(i,:), shuttle_data, Qs)));
                new_violation = sum(max(0, c_new));
                if new_violation < old_violation || new_fval < fitness(i)
                    update_individual = true;
                end
            end
            
            if update_individual
                fitness(i) = new_fval;
                is_feasible(i) = new_feasible;
            end
        end
        
        % ========== 记录历史 ==========
        all_params_history = [all_params_history; population];
        all_fvals_history = [all_fvals_history; fitness];
        
        % ========== 显示进度 ==========
        if mod(iter, 20) == 0 || iter <= 5
            best_feasible_idx = find_best_feasible(fitness, is_feasible);
            fprintf('  迭代 %d/%d, R2=%.2f, ST=%.2f, 可行解比例: %.1f%%, 最佳适应度: %.4f\n', ...
                iter, max_iterations, R2, ST, feasibility_history(iter)*100, fitness(best_feasible_idx));
        end
        
        % ========== 动态调整参数 ==========
        if mod(iter, 10) == 0
            % 调整安全阈值
            if feasible_count < 0.2 * n_pop
                ST = min(0.9, ST + 0.05);  % 增加安全阈值
            else
                ST = max(0.3, ST - 0.05);  % 降低安全阈值
            end
        end
    end
    
    % ========== 找到最佳参数 ==========
    [best_fval, best_idx] = min_fitness_with_feasibility(fitness, is_feasible);
    best_params = population(best_idx, :);
    
    % 验证最佳参数
    [c_best, ~] = parameter_constraints_with_all_metrics_3504(best_params, shuttle_data, Qs);
    best_feasible = ~any(c_best > 0);
    
    % 保存结果
    results = struct();
    results.best_fval = best_fval;
    results.all_params = all_params_history;
    results.all_fvals = all_fvals_history;
    results.iterations = max_iterations;
    results.feasibility_history = feasibility_history;
    results.is_feasible = best_feasible;
    
    fprintf('增强的约束增强麻雀搜索算法完成！最佳目标函数值: %f, 可行: %s\n', ...
        best_fval, string(best_feasible));
end

% 针对3504组数据的多起始点fmincon优化（使用所有起点）
function [best_params, best_fval, info] = run_constrained_multi_start_fmincon_for_method_3504(...
        starting_points, shuttle_data, Qs, lb, ub, fmincon_options, constraint_weight)
    
    % 使用所有起点，不限制数量
    n_points = size(starting_points, 1);
    if n_points < 100
        fprintf('  警告：只有%d个起点，少于100个\n', n_points);
    end
    
    fprintf('  使用 %d 个起始点运行约束增强fmincon...\n', n_points);
    
    % 定义约束增强的目标函数
    fmincon_constrained_objective = @(x) unified_constrained_objective_3504(x, shuttle_data, Qs, constraint_weight);
    
    % 约束函数
    nonlcon = @(x) parameter_constraints_with_all_metrics_3504(x, shuttle_data, Qs);
    
    % 预分配结果存储
    all_params_cell = cell(n_points, 1);
    all_fvals = zeros(n_points, 1);
    all_is_feasible = false(n_points, 1);
    all_exitflags = zeros(n_points, 1);
    all_success = false(n_points, 1);
    
 
    for i = 1:n_points
        try
            [params, fval, exitflag, output] = ...
                fmincon(fmincon_constrained_objective, ...
                        starting_points(i, :), [], [], [], [], lb, ub, [], fmincon_options);
            
            % 因为去掉了 nonlcon，所有的解我们都先当作“软可行解”
            is_feasible = true; 
            
            all_params_cell{i} = params;
            all_fvals(i) = fval;
            all_is_feasible(i) = is_feasible;
            all_exitflags(i) = exitflag;
            all_success(i) = true;
            
        catch ME
            fval = 100; % 平滑惩罚
            all_params_cell{i} = starting_points(i, :);
            all_fvals(i) = fval;
            all_is_feasible(i) = false;
            all_exitflags(i) = -1;
            all_success(i) = false;
        end
        
        if mod(i, 10) == 0
            fprintf('    完成 %d/%d 个起始点\n', i, n_points);
        end
    end
    
    % 将cell数组转换为矩阵
    all_params_matrix = zeros(n_points, length(lb));
    for i = 1:n_points
        if ~isempty(all_params_cell{i})
            all_params_matrix(i, :) = all_params_cell{i};
        else
            all_params_matrix(i, :) = starting_points(i, :);
        end
    end
    
    % 选择目标函数值最小的解
    [best_fval_local, best_idx] = min(all_fvals);
    best_params_local = all_params_matrix(best_idx, :);
    
    % ===================== 【修改点2】：精英保留机制 (Elitism) =====================
    % 评估传入的 100 个历史起点在当前 KGE 目标函数下的真实得分
    starting_fvals = zeros(n_points, 1);
    for k = 1:n_points
        starting_fvals(k) = unified_constrained_objective_3504(starting_points(k, :), shuttle_data, Qs, constraint_weight);
    end
    [best_starting_fval, best_start_idx] = min(starting_fvals);

    % 如果 fmincon 陷入了劣质针眼导致性能不如传入的最佳起点，强行回滚！
    if best_fval_local > best_starting_fval
        fprintf('  [精英保留]: fmincon 微调后性能退化 (%.4f > %.4f)，回滚至单一算法原始最优解。\n', best_fval_local, best_starting_fval);
        best_fval_local = best_starting_fval;
        best_params_local = starting_points(best_start_idx, :);
    end
    % ==============================================================================
    
    % 统计信息
    success_count = sum(all_success);
    info = struct(...
        'num_starting_points', n_points, ...
        'success_count', success_count, ...
        'all_fvals', all_fvals, ...
        'min_fval', min(all_fvals), ...
        'all_params', all_params_matrix);
    
    best_params = best_params_local;
    best_fval = best_fval_local;
end


% ===================== 新增：质量过滤函数 =====================
function [filtered_params, filtered_fvals] = filter_high_quality_params(params, fvals, lb, ub, shuttle_data, Qs, constraint_weight, min_quality_ratio)
    % 过滤高质量参数
    % min_quality_ratio: 保留前多少比例的参数（0-1）
    
    if nargin < 8
        min_quality_ratio = 0.50 ;  % 默认保留前50%
    end
    
    % 计算质量分数（目标函数值越小越好）
    normalized_fvals = (fvals - min(fvals)) / (max(fvals) - min(fvals) + eps);
    quality_scores = 1 - normalized_fvals;  % 转换为质量分数（越高越好）
    
    % 计算分位数阈值
    threshold = prctile(quality_scores, 100 * (1 - min_quality_ratio));
    
    % 选择高质量参数
    high_quality_idx = quality_scores >= threshold;
    filtered_params = params(high_quality_idx, :);
    filtered_fvals = fvals(high_quality_idx);
end

% ===================== 通用辅助函数 =====================

% 计算纳什-萨特克利夫效率系数
function nse = calculate_NSE(sim, obs)
    valid_idx = ~isnan(sim) & ~isnan(obs);
    sim = sim(valid_idx);
    obs = obs(valid_idx);
    
    if isempty(sim) || isempty(obs) || length(sim) < 2
        nse = -Inf;
        return;
    end
    
    numerator = sum((sim - obs).^2);
    denominator = sum((obs - mean(obs)).^2);
    
    if denominator == 0
        nse = -Inf;
    else
        nse = 1 - (numerator / denominator);
    end
end

% ========== 新增：拉丁超立方采样函数 ==========
function samples = lhsdesign(n_samples, n_dim)
    % 拉丁超立方采样
    samples = zeros(n_samples, n_dim);
    for i = 1:n_dim
        samples(:, i) = (randperm(n_samples)' - rand(n_samples, 1)) / n_samples;
    end
end

% ========== 新增：Moore-Penrose伪逆函数 ==========
function A_plus = pinv(A)
    % 简化的伪逆计算（避免使用MATLAB的pinv函数）
    [U, S, V] = svd(A);
    tolerance = max(size(A)) * eps(norm(S, inf));
    s = diag(S);
    s(s > tolerance) = 1 ./ s(s > tolerance);
    S_plus = diag(s);
    A_plus = V * S_plus * U';
end

% ========== 增强的排序函数（考虑可行性和多样性） ==========
function [sorted_fitness, sorted_idx] = sort_fitness_with_feasibility_diversity(fitness, feasibility, positions)
    n = length(fitness);
    compound_fitness = zeros(n, 1);
    
    for i = 1:n
        if feasibility(i)
            % 可行解有优势
            compound_fitness(i) = fitness(i) - 1e6;
        else
            % 不可行解，但考虑多样性
            diversity_score = calculate_diversity(positions(i,:), positions);
            compound_fitness(i) = fitness(i) - 0.1 * diversity_score;
        end
    end
    
    [~, sorted_idx] = sort(compound_fitness);
    sorted_fitness = fitness(sorted_idx);
end

function diversity = calculate_diversity(individual, population)
    % 计算个体与种群的距离（多样性）
    n = size(population, 1);
    distances = zeros(n, 1);
    
    for i = 1:n
        distances(i) = norm(individual - population(i,:));
    end
    
    diversity = mean(distances);
end

% 考虑可行性的适应度排序
function [sorted_fitness, sorted_idx] = sort_fitness_with_feasibility(fitness, is_feasible)
    n = length(fitness);
    compound_fitness = zeros(n, 1);
    for i = 1:n
        if is_feasible(i)
            compound_fitness(i) = fitness(i) - 1e8;  % 可行解优先
        else
            compound_fitness(i) = fitness(i);
        end
    end
    
    [~, sorted_idx] = sort(compound_fitness);
    sorted_fitness = fitness(sorted_idx);
end

% 考虑可行性的最小适应度查找
function [min_fitness, min_idx] = min_fitness_with_feasibility(fitness, is_feasible)
    feasible_indices = find(is_feasible);
    
    if ~isempty(feasible_indices)
        [min_fitness, min_idx_in_feasible] = min(fitness(feasible_indices));
        min_idx = feasible_indices(min_idx_in_feasible);
    else
        [min_fitness, min_idx] = min(fitness);
    end
end

% 找到最佳可行个体
function best_idx = find_best_feasible(fitness, feasibility)
    feasible_indices = find(feasibility);
    
    if ~isempty(feasible_indices)
        [~, best_idx_in_feasible] = min(fitness(feasible_indices));
        best_idx = feasible_indices(best_idx_in_feasible);
    else
        [~, best_idx] = min(fitness);
    end
end
% ===================== 参数性能约束评估函数 =====================
function [c, ceq] = parameter_constraints_with_all_metrics_3504(params, shuttle_data, Qs)
    % 非线性不等式约束 c <= 0 （满足约束时 c 必须小于等于0）
    % 非线性等式约束 ceq = 0
    ceq = [];
    c = zeros(1, 4);
    
    try
        % 运行模型获取各项统计指标
        [is_valid, stats] = validate_parameters_with_all_metrics_3504(params, shuttle_data, Qs);
        
        if is_valid && stats.ET_samples >= 50
            % 设置硬性约束 (c <= 0 表示合规)
            c(1) = 0.9 - stats.slope_ET;  % 要求斜率 slope_ET >= 0.9
            c(2) = stats.slope_ET - 1.1;  % 要求斜率 slope_ET <= 1.1
            c(3) = 0.80 - stats.R2_ET;     % 要求 R2_ET >= 0.80
            c(4) = 0.80 - stats.NSE_ET;    % 要求 NSE_ET >= 0.80
            
            % 如果出现计算错误导致NaN，直接赋大正数表示严重违规
            if any(isnan(c))
                c = [100, 100, 100, 100];
            end
        else
            % 模型失效或者有效样本过少，严重违规
            c = [100, 100, 100, 100];
        end
    catch
        % 代码崩溃时，严重违规
        c = [100, 100, 100, 100];
    end
end