% =========================================================================
% 终极自动化后处理与统计脚本 v5.1 (仅提取: BO-Fmincon + 策略 B)
% =========================================================================
clearvars -except shuttle_FlagBadData history Qs shuttle_original all_method_params shuttle_calibration shuttle_carlibration;
clc;

% ====== 数据加载与基础设置 ======
if ~exist('history', 'var')
    if exist('optimization_results_3504.mat', 'file')
        load('optimization_results_3504.mat', 'history', 'shuttle_FlagBadData', 'Qs');
        fprintf('已从 optimization_results_3504.mat 加载历史数据和气象数据。\n');
    else
        error('未找到 history 变量且 optimization_results_3504.mat 文件不存在。请先运行主优化程序。');
    end
end
if ~exist('shuttle_FlagBadData', 'var'), error('工作区中未找到 shuttle_FlagBadData 数据，请检查 mat 文件内容。'); end
if ~exist('Qs', 'var'), Qs = 0.46; end

% 定义算法与策略
method_list = {'bayes', 'ga', 'pso', 'ssa'};
algorithms = {}; 
for m_idx = 1:length(method_list)
    method = method_list{m_idx};
    orig_fvals = history.all_methods.(method).all_fvals;
    orig_params = history.all_methods.(method).all_params;
    [~, sort_idx] = sort(orig_fvals);
    algorithms{end+1} = struct('name', method, 'top10', orig_params(sort_idx(1:min(10, end)), :));
    method_metrics = history.all_method_metrics.(method);
    algorithms{end+1} = struct('name', [method '-fmincon'], 'top10', method_metrics.top10_params);
end
strategies = {'A', 'B', 'C'}; 
years = [2020, 2022];

% ============= 新增：专门用来收集 BO-Fmincon+策略B 结果的容器 =============
bo_result_table = {};

% ====== 核心循环 ======
for a_idx = 1:length(algorithms)
    alg_name = algorithms{a_idx}.name;
    top10_params = algorithms{a_idx}.top10;
    
    for s_idx = 1:3
        strategy_name = strategies{s_idx};
        
        % 【拦截条件】只处理 BO-Fmincon (即 bayes-fmincon) 且 策略为 B 的数据
        % 如果您的变量名是 'BO-fmincon'，请把下方 'bayes-fmincon' 改成 'BO-fmincon'
        if ~(strcmp(alg_name, 'bayes-fmincon') && strcmp(strategy_name, 'B'))
            continue; 
        end
        
        fprintf('正在提取 BO-Fmincon + 策略 B 的生育期数据...\n');
        
        if s_idx == 1
            param_A = make_para_matrix(top10_params(1, :));
            hh_matrix_raw = halfhour_shuttleworth_fillgap(shuttle_FlagBadData, param_A, Qs);
        elseif s_idx == 2
            param_B = make_para_matrix(mean(top10_params, 1));
            hh_matrix_raw = halfhour_shuttleworth_fillgap(shuttle_FlagBadData, param_B, Qs);
        elseif s_idx == 3
            T_ens = zeros(size(shuttle_FlagBadData, 1), 1); E_ens = zeros(size(shuttle_FlagBadData, 1), 1); ET_ens = zeros(size(shuttle_FlagBadData, 1), 1);
            for k = 1:size(top10_params, 1)
                param_k = make_para_matrix(top10_params(k, :));
                temp_hh = halfhour_shuttleworth_fillgap(shuttle_FlagBadData, param_k, Qs);
                T_ens = T_ens + temp_hh(:,21); E_ens = E_ens + temp_hh(:,22); ET_ens = ET_ens + temp_hh(:,23);
            end
            T_ens = T_ens / size(top10_params, 1); E_ens = E_ens / size(top10_params, 1); ET_ens = ET_ens / size(top10_params, 1);
            hh_matrix_raw = temp_hh; 
            hh_matrix_raw(:,21) = T_ens; hh_matrix_raw(:,22) = E_ens; hh_matrix_raw(:,23) = ET_ens;
            hh_matrix_raw(:,24) = E_ens ./ ET_ens; 
        end
        
        for y = 1:length(years)
            yr = years(y);
            idx_yr = (hh_matrix_raw(:,1) == yr);
            hh_yr_raw = hh_matrix_raw(idx_yr, :);
            if isempty(hh_yr_raw), continue; end
            
            [daily_yr, ~] = shuttlesum(hh_yr_raw, ''); 
            
            % ========== 核心：基于东区农学记录的 4 个精细生育期提取 ==========
            DOY = daily_yr(:, 4);
            if yr == 2020
                idx_emg = DOY >= 184 & DOY < 213;  % 1. 苗期 
                idx_jnt = DOY >= 213 & DOY < 231;  % 2. 拔节期 
                idx_tas = DOY >= 231 & DOY < 259;  % 3. 抽穗期 
                idx_fil = DOY >= 259 & DOY <= 269; % 4. 灌浆成熟期 
            else
                idx_emg = DOY >= 184 & DOY < 207;  % 1. 苗期 
                idx_jnt = DOY >= 207 & DOY < 228;  % 2. 拔节期 
                idx_tas = DOY >= 228 & DOY < 257;  % 3. 抽穗期 
                idx_fil = DOY >= 257 & DOY <= 268; % 4. 灌浆成熟期 
            end
            
            % 封装阶段名称与索引（方便后续循环输出）
            stage_idxs = {idx_emg, idx_jnt, idx_tas, idx_fil};
            stage_names = {'苗期', '拔节期', '抽穗期', '灌浆成熟期'};
            
            % 提取指定 6 列数据：生育期，年份，天数，实测ET均值，模拟ET均值，阶段内RMSE
            % 注：实测ET列号=18，模拟ET列号=22
            for st = 1:4
                cur_idx = stage_idxs{st};
                days = sum(cur_idx);
                if days > 0
                    obs_ET_stage = daily_yr(cur_idx, 18);
                    sim_ET_stage = daily_yr(cur_idx, 22);
                    
                    mean_obs = mean(obs_ET_stage, 'omitnan');
                    mean_sim = mean(sim_ET_stage, 'omitnan');
                    rmse_stage = sqrt(mean((obs_ET_stage - sim_ET_stage).^2, 'omitnan'));
                else
                    mean_obs = NaN; mean_sim = NaN; rmse_stage = NaN;
                end
                
                % 将本次结果存入容器
                bo_result_table = [bo_result_table; {stage_names{st}, yr, days, mean_obs, mean_sim, rmse_stage}];
            end
        end
    end
end

% ====== 最终：生成并打印整洁的表格 ======
if ~isempty(bo_result_table)
    % 转换成 MATLAB 原生 Table，方便直接在控制台预览
    T_bo = cell2table(bo_result_table, ...
        'VariableNames', {'生育期', '年份', '阶段内天数', 'ET实测均值', 'ET模拟均值(BO-Fmincon)', '阶段内RMSE'});
    
    fprintf('\n=======================================================\n');
    fprintf('  BO-Fmincon + 策略 B 生育期评价完整结果 \n');
    fprintf('=======================================================\n');
    disp(T_bo);
    
    % 如果想直接保存成 Excel 文件使用，请取消下面这行的注释
    % writetable(T_bo, 'BO-Fmincon_策略B_生育期结果.xlsx');
else
    fprintf('没有提取到 BO-Fmincon 或策略 B 的数据，请检查 alg_name 是否写对。\n');
end

% ====== 必须保证原脚本需要的辅助函数存在 ======
function p_matrix = make_para_matrix(p_vec)
    p_matrix = zeros(4,6); p_matrix(3,:) = p_vec; 
end