% =========================================================================
% Generate_Strategy_Table.m
% 提取 8 种算法（原始 + fmincon 混合）的三种策略（A、B、C）参数及 SGPI
% =========================================================================
clear; clc;

% 1. 加载优化结果
if ~exist('optimization_results_3504.mat', 'file')
    error('未找到 optimization_results_3504.mat，请确保文件在当前目录。');
end
load('optimization_results_3504.mat', 'history', 'shuttle_carlibration', 'Qs');

% 若 shuttle_carlibration 未保存，则尝试从 history 中提取（以防万一）
if ~exist('shuttle_carlibration', 'var')
    error('未找到 shuttle_carlibration，请检查 mat 文件内容。');
end

% 2. 定义 8 种算法名称
algorithms = {'bayes', 'ga', 'pso', 'ssa', 'bayes_f', 'ga_f', 'pso_f', 'ssa_f'};

% 预存结果的结构体数组
Results = struct('Algorithm', {}, 'Strategy', {}, 'b1', {}, 'b2', {}, 'b3', {}, ...
                 'a1', {}, 'g01', {}, 'LightExtCoef', {}, 'SGPI', {});

% 3. 遍历每种算法
for alg_idx = 1:length(algorithms)
    alg = algorithms{alg_idx};
    is_hybrid = contains(alg, '_f');
    base_alg = strrep(alg, '_f', '');
    
    fprintf('处理算法: %s\n', alg);
    
    % ----- 提取该算法的候选参数集（用于策略 A, B, C）-----
    if ~is_hybrid
        % 原始算法：直接从 history.all_methods 获取
        if ~isfield(history, 'all_methods') || ~isfield(history.all_methods, base_alg)
            warning('未找到算法 %s 的历史数据，跳过。', base_alg);
            continue;
        end
        all_p = history.all_methods.(base_alg).all_params;
        all_f = history.all_methods.(base_alg).all_fvals;
        % 按 fval 升序排序
        [sorted_f, sort_idx] = sort(all_f);
        sorted_p = all_p(sort_idx, :);
        % 取前 10 个（若不足则全取）
        n_top = min(10, size(sorted_p, 1));
        top10_params = sorted_p(1:n_top, :);
        top10_fvals = sorted_f(1:n_top);
    else
        % 混合算法：从 history.all_method_metrics 中提取 is_fmincon_optimized == true 的参数
        if ~isfield(history, 'all_method_metrics') || ~isfield(history.all_method_metrics, base_alg)
            warning('未找到算法 %s 的 metrics，跳过。', base_alg);
            continue;
        end
        metrics = history.all_method_metrics.(base_alg);
        if ~isfield(metrics, 'top10_params') || ~isfield(metrics, 'top10_fvals') || ~isfield(metrics, 'is_fmincon_optimized')
            warning('算法 %s 的 metrics 缺少所需字段，跳过。', base_alg);
            continue;
        end
        % 提取所有 fmincon 优化的参数（标记为 true）
        all_hybrid_p = metrics.top10_params(metrics.is_fmincon_optimized, :);
        all_hybrid_f = metrics.top10_fvals(metrics.is_fmincon_optimized);
        if isempty(all_hybrid_p)
            warning('算法 %s 没有 fmincon 优化结果，跳过。', base_alg);
            continue;
        end
        % 按 fval 排序
        [sorted_f, sort_idx] = sort(all_hybrid_f);
        sorted_p = all_hybrid_p(sort_idx, :);
        % 取前 10 个
        n_top = min(10, size(sorted_p, 1));
        top10_params = sorted_p(1:n_top, :);
        top10_fvals = sorted_f(1:n_top);
    end
    
    % 如果有效参数不足 1，则跳过
    if size(top10_params, 1) < 1
        warning('算法 %s 没有可用参数，跳过。', alg);
        continue;
    end
    
    % ----- 策略 A：最优单组 -----
    params_A = top10_params(1, :);
    fval_A = top10_fvals(1);
    
    % ----- 策略 B：参数平均（前 n_top 组）-----
    params_B = mean(top10_params, 1);
    % 计算策略 B 的 SGPI（重新评估均值参数的目标值）
    fval_B = calc_SGPI(params_B, shuttle_carlibration, Qs);
    
    % ----- 策略 C：集合平均（前 n_top 组回代，取 ET 平均）-----
    % 初始化累加 ET
    ET_ensemble = zeros(size(shuttle_carlibration, 1), 1);
    for k = 1:size(top10_params, 1)
        p = top10_params(k, :);
        % 调用 SWH 模型（需确保 SWH_halfhour 在路径中）
        [shuttle_temp] = SWH_halfhour(shuttle_carlibration, p(1), p(2), p(3), p(4), p(5), p(6), Qs);
        ET_sim = shuttle_temp(:, 23);
        ET_ensemble = ET_ensemble + ET_sim;
    end
    ET_ensemble = ET_ensemble / size(top10_params, 1);
    % 提取实测 ET
    ET_obs = shuttle_carlibration(:, 17);
    % 计算策略 C 的 SGPI
    fval_C = calc_SGPI_from_vectors(ET_ensemble, ET_obs);
    
    % ----- 存储结果 -----
    % 策略 A
    Results(end+1) = struct('Algorithm', alg, 'Strategy', 'A', ...
        'b1', params_A(1), 'b2', params_A(2), 'b3', params_A(3), ...
        'a1', params_A(4), 'g01', params_A(5), 'LightExtCoef', params_A(6), ...
        'SGPI', fval_A);
    % 策略 B
    Results(end+1) = struct('Algorithm', alg, 'Strategy', 'B', ...
        'b1', params_B(1), 'b2', params_B(2), 'b3', params_B(3), ...
        'a1', params_B(4), 'g01', params_B(5), 'LightExtCoef', params_B(6), ...
        'SGPI', fval_B);
    % 策略 C
    Results(end+1) = struct('Algorithm', alg, 'Strategy', 'C', ...
        'b1', NaN, 'b2', NaN, 'b3', NaN, ...  % 策略 C 没有单一参数，留空
        'a1', NaN, 'g01', NaN, 'LightExtCoef', NaN, ...
        'SGPI', fval_C);
end

% 4. 将结构体转为表格并写入 Excel
T = struct2table(Results);
writetable(T, 'Strategy_Table.xlsx');
fprintf('\n? 成功生成表格 Strategy_Table.xlsx，包含 %d 行数据。\n', height(T));

% ========== 辅助函数 ==========

function sgpi = calc_SGPI(params, shuttle_data, Qs)
    % 计算单组参数的 SGPI（与目标函数一致）
    try
        [shuttle] = SWH_halfhour(shuttle_data, params(1), params(2), params(3), ...
                                  params(4), params(5), params(6), Qs);
        ET_sim = shuttle(:,23);
        ET_obs = shuttle(:,17);
        sgpi = compute_SGPI(ET_sim, ET_obs);
    catch
        sgpi = NaN;
    end
end

function sgpi = calc_SGPI_from_vectors(ET_sim, ET_obs)
    % 直接传入模拟和观测向量计算 SGPI
    sgpi = compute_SGPI(ET_sim, ET_obs);
end

function sgpi = compute_SGPI(sim, obs)
    % 计算 SGPI = (1-R2) + (1-NSE) + RMSE/mean_obs + MAE/mean_obs
    valid = (obs ~= -99999) & ~isnan(obs) & ~isnan(sim) & ~isinf(sim);
    sim = sim(valid);
    obs = obs(valid);
    if length(sim) < 10
        sgpi = NaN;
        return;
    end
    % R2
    cc = corrcoef(sim, obs);
    r2 = cc(1,2)^2;
    % NSE
    nse = 1 - sum((sim - obs).^2) / sum((obs - mean(obs)).^2);
    % RMSE, MAE
    rmse = sqrt(mean((sim - obs).^2));
    mae = mean(abs(sim - obs));
    mean_obs = mean(obs);
    if mean_obs == 0, mean_obs = eps; end
    sgpi = (1 - r2) + (1 - nse) + (rmse / mean_obs) + (mae / mean_obs);
end