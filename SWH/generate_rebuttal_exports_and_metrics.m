% =========================================================================
% 终极自动化后处理与统计脚本 v5.1
% 新增: Bias (物理偏差方向), KGE, 以及截距 (Intercept)
% =========================================================================
clearvars -except shuttle_FlagBadData history Qs shuttle_original all_method_params shuttle_calibration shuttle_carlibration;
clc;

fprintf('=======================================================\n');
fprintf(' 启动批处理 v5.1: 新增 Bias, KGE, 截距 及多尺度组内 GPI 严格排名\n');
fprintf('=======================================================\n\n');

% ====== 【新增】若工作区缺少 history，则从 .mat 文件加载 ======
if ~exist('history', 'var')
    if exist('optimization_results_3504.mat', 'file')
        load('optimization_results_3504.mat', 'history', 'shuttle_FlagBadData', 'Qs');
        fprintf('已从 optimization_results_3504.mat 加载历史数据和气象数据。\n');
    else
        error('未找到 history 变量且 optimization_results_3504.mat 文件不存在。请先运行主优化程序。');
    end
end

if ~exist('shuttle_FlagBadData', 'var')
    error('工作区中未找到 shuttle_FlagBadData 数据，请检查 mat 文件内容。');
end
if ~exist('Qs', 'var'), Qs = 0.46; end

if ~exist('Export_2020', 'dir'), mkdir('Export_2020'); end
if ~exist('Export_2022', 'dir'), mkdir('Export_2022'); end

raw_results = []; 
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

% ---------------- 阶段一：高精数据生成与单位换算 ----------------
for a_idx = 1:length(algorithms)
    alg_name = algorithms{a_idx}.name;
    top10_params = algorithms{a_idx}.top10;
    fprintf('正在处理算法并计算新增指标(Bias, KGE, Intercept): %-15s (%d/8)...\n', alg_name, a_idx);
    
    for s_idx = 1:3
        strategy_name = strategies{s_idx};
        
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
            
            hh_yr_export = hh_yr_raw;
            cols_to_convert = [17, 21, 22, 23];
            for c = cols_to_convert
                valid = hh_yr_export(:, c) ~= -99999;
                hh_yr_export(valid, c) = hh_yr_export(valid, c) * 3.6;
            end
            
            % 半小时尺度指标 (新增 Intercept)
            [K_hh, Intercept_hh, R2_hh, RMSE_hh, NSE_hh, MAE_hh, Bias_hh, KGE_hh] = calc_metrics(hh_yr_export(:,23), hh_yr_export(:,17));
            res_hh = struct('Year', yr, 'Scale', 'Half-Hourly', 'Algorithm', alg_name, 'Strategy', strategy_name, ...
                            'K', K_hh, 'Intercept', Intercept_hh, 'R2', R2_hh, 'RMSE', RMSE_hh, 'NSE', NSE_hh, 'MAE', MAE_hh, 'Bias', Bias_hh, 'KGE', KGE_hh, 'BaseAlgo', strrep(alg_name, '-fmincon', ''));
            raw_results = [raw_results; res_hh];
            
            % 日尺度指标 (新增 Intercept)
            [K_d, Intercept_d, R2_d, RMSE_d, NSE_d, MAE_d, Bias_d, KGE_d] = calc_metrics(daily_yr(:,22), daily_yr(:,18));
            % ========== 新增：基于东区农学记录的 4 个精细生育期验证 ==========
            DOY = daily_yr(:, 4);
            if yr == 2020
                % 2020 (闰年) 节点: 出苗184(7/2), 拔节213(7/31), 抽穗231(8/18), 乳熟259(9/15), 测产269(9/25)
                idx_emg = DOY >= 184 & DOY < 213;  % 1. 苗期 (出苗至拔节前)
                idx_jnt = DOY >= 213 & DOY < 231;  % 2. 拔节期 (拔节至抽穗前)
                idx_tas = DOY >= 231 & DOY < 259;  % 3. 抽穗期 (抽穗至乳熟前，耗水绝对高峰)
                idx_fil = DOY >= 259 & DOY <= 269; % 4. 灌浆成熟期 (乳熟至测产，严防裸地数据)
            else
                % 2022 (平年-东区) 节点: 出苗184(7/3), 拔节207(7/26), 抽穗228(8/16), 乳熟257(9/14), 成熟268(9/25)
                idx_emg = DOY >= 184 & DOY < 207;  % 1. 苗期 
                idx_jnt = DOY >= 207 & DOY < 228;  % 2. 拔节期 
                idx_tas = DOY >= 228 & DOY < 257;  % 3. 抽穗期 
                idx_fil = DOY >= 257 & DOY <= 268; % 4. 灌浆成熟期 
            end
            
            % 提取各阶段的 4 大核心指标
            [~, ~, R2_emg, RMSE_emg, NSE_emg, ~, ~, KGE_emg] = calc_metrics(daily_yr(idx_emg,22), daily_yr(idx_emg,18));
            [~, ~, R2_jnt, RMSE_jnt, NSE_jnt, ~, ~, KGE_jnt] = calc_metrics(daily_yr(idx_jnt,22), daily_yr(idx_jnt,18));
            [~, ~, R2_tas, RMSE_tas, NSE_tas, ~, ~, KGE_tas] = calc_metrics(daily_yr(idx_tas,22), daily_yr(idx_tas,18));
            [~, ~, R2_fil, RMSE_fil, NSE_fil, ~, ~, KGE_fil] = calc_metrics(daily_yr(idx_fil,22), daily_yr(idx_fil,18));
            
            % 打印基于东区对齐的精细化报表
            fprintf('  [%d年 农学精细分期验证结果 (R? | RMSE | NSE | KGE)]\n', yr);
            fprintf('   -> 1. 苗期     (出苗-拔节): R?=%.3f, RMSE=%.3f mm/d, NSE=%.3f, KGE=%.3f\n', R2_emg, RMSE_emg, NSE_emg, KGE_emg);
            fprintf('   -> 2. 拔节期   (拔节-抽穗): R?=%.3f, RMSE=%.3f mm/d, NSE=%.3f, KGE=%.3f\n', R2_jnt, RMSE_jnt, NSE_jnt, KGE_jnt);
            fprintf('   -> 3. 抽穗期   (抽穗-乳熟): R?=%.3f, RMSE=%.3f mm/d, NSE=%.3f, KGE=%.3f\n', R2_tas, RMSE_tas, NSE_tas, KGE_tas);
            fprintf('   -> 4. 灌浆成熟 (乳熟-测产): R?=%.3f, RMSE=%.3f mm/d, NSE=%.3f, KGE=%.3f\n', R2_fil, RMSE_fil, NSE_fil, KGE_fil);
            fprintf('  --------------------------------------------------------------\n');
            % ==============================================================================
            res_d = struct('Year', yr, 'Scale', 'Daily', 'Algorithm', alg_name, 'Strategy', strategy_name, ...
                           'K', K_d, 'Intercept', Intercept_d, 'R2', R2_d, 'RMSE', RMSE_d, 'NSE', NSE_d, 'MAE', MAE_d, 'Bias', Bias_d, 'KGE', KGE_d, 'BaseAlgo', strrep(alg_name, '-fmincon', ''));
            raw_results = [raw_results; res_d];
            
            export_data_files_standard(hh_yr_export, daily_yr, yr, alg_name, strategy_name);
        end
    end
end

% ---------------- 阶段二：组内归一化、GPI及各单项辅助排名分析 ----------------
fprintf('\n正在执行组内标准化、多维度辅助排名及核心指标提升度统计...\n');

% 严格按照用户要求的顺序更新表头（新增 Intercept）
metrics_summary = {'Year', 'Scale', 'Algorithm', 'Strategy', 'K(Slope)', 'Intercept', 'R2', 'RMSE', 'NSE', 'MAE', 'Bias', 'KGE', ...
                   'Rank_R2', 'Rank_NSE', 'Rank_RMSE', 'Rank_MAE', 'Rank_absBias', 'Rank_KGE', ...
                   'GPI_Score', 'Rank_GPI', ...
                   'R2_Change(‰)', 'NSE_Change(‰)', 'RMSE_Change(‰)', 'MAE_Change(‰)'};
scales = {'Half-Hourly', 'Daily'};

for y = 1:length(years)
    for s = 1:length(scales)
        yr = years(y);
        sc = scales{s};
        
        idx = find([raw_results.Year] == yr & strcmp({raw_results.Scale}, sc));
        if isempty(idx), continue; end
        group_data = raw_results(idx);
        
        % 提取组内数组
        K_arr = [group_data.K]; R2_arr = [group_data.R2]; RMSE_arr = [group_data.RMSE];
        NSE_arr = [group_data.NSE]; MAE_arr = [group_data.MAE]; KGE_arr = [group_data.KGE]; Bias_arr = [group_data.Bias];
        
        % ========== 1. 计算单项指标的独立严格排名 (补全了 R2 和 MAE) ==========
        [~, sort_r2] = sort(R2_arr, 'descend'); rank_R2 = zeros(size(R2_arr)); rank_R2(sort_r2) = 1:length(R2_arr);
        [~, sort_nse] = sort(NSE_arr, 'descend'); rank_NSE = zeros(size(NSE_arr)); rank_NSE(sort_nse) = 1:length(NSE_arr);
        [~, sort_rmse] = sort(RMSE_arr, 'ascend'); rank_RMSE = zeros(size(RMSE_arr)); rank_RMSE(sort_rmse) = 1:length(RMSE_arr);
        [~, sort_mae] = sort(MAE_arr, 'ascend'); rank_MAE = zeros(size(MAE_arr)); rank_MAE(sort_mae) = 1:length(MAE_arr);
        [~, sort_bias] = sort(abs(Bias_arr), 'ascend'); rank_absBias = zeros(size(Bias_arr)); rank_absBias(sort_bias) = 1:length(Bias_arr);
        [~, sort_kge] = sort(KGE_arr, 'descend'); rank_KGE = zeros(size(KGE_arr)); rank_KGE(sort_kge) = 1:length(KGE_arr);
        
        % ========== 2. 依据 Jia et al. (2022) 的中位数偏差法计算 GPI ==========
        scale_minmax = @(x) (x - min(x)) ./ (max(x) - min(x) + eps);
        y_R2 = scale_minmax(R2_arr);
        y_NSE = scale_minmax(NSE_arr);
        y_RMSE = scale_minmax(RMSE_arr);
        y_MAE = scale_minmax(MAE_arr);
        
        g_R2 = median(y_R2);
        g_NSE = median(y_NSE);
        g_RMSE = median(y_RMSE);
        g_MAE = median(y_MAE);
        
        % alpha = -1 (正向指标), alpha = 1 (负向指标)
        GPI_arr = -1 * (g_R2 - y_R2) + -1 * (g_NSE - y_NSE) + 1 * (g_RMSE - y_RMSE) + 1 * (g_MAE - y_MAE);
                   
        [~, sort_gpi] = sort(GPI_arr, 'descend'); 
        rank_GPI = zeros(size(GPI_arr)); 
        rank_GPI(sort_gpi) = 1:length(GPI_arr);
        
        % ========== 3. 收集并按照严格顺序写入报表 ==========
        for i = 1:length(group_data)
            alg = group_data(i).Algorithm;
            r2_chg = NaN; rmse_chg = NaN; nse_chg = NaN; mae_chg = NaN;
            
            % 计算千分比变化 (相对于基础算法)
            if contains(alg, '-fmincon')
                base_alg = group_data(i).BaseAlgo;
                strat = group_data(i).Strategy;
                base_idx = find(strcmp({group_data.Algorithm}, base_alg) & strcmp({group_data.Strategy}, strat));
                if ~isempty(base_idx)
                    base_res = group_data(base_idx(1));
                    r2_chg = (group_data(i).R2 - base_res.R2) / abs(base_res.R2) * 1000;
                    nse_chg = (group_data(i).NSE - base_res.NSE) / abs(base_res.NSE) * 1000;
                    rmse_chg = (group_data(i).RMSE - base_res.RMSE) / abs(base_res.RMSE) * 1000;
                    mae_chg = (group_data(i).MAE - base_res.MAE) / abs(base_res.MAE) * 1000;
                end
            end
            
            % 严格对齐列顺序（新增 Intercept）
            metrics_summary(end+1, :) = {yr, sc, alg, group_data(i).Strategy, ...
                                         group_data(i).K, group_data(i).Intercept, group_data(i).R2, group_data(i).RMSE, ...
                                         group_data(i).NSE, group_data(i).MAE, group_data(i).Bias, group_data(i).KGE, ...
                                         rank_R2(i), rank_NSE(i), rank_RMSE(i), rank_MAE(i), rank_absBias(i), rank_KGE(i), ...
                                         GPI_arr(i), rank_GPI(i), ...
                                         r2_chg, nse_chg, rmse_chg, mae_chg};
        end
    end
end

try
    xlswrite('Summary_Metrics_Full_Report_with_Intercept.xlsx', metrics_summary, 'Sheet1', 'A1');
    fprintf('=> 成功生成严格排序排版的报表: Summary_Metrics_Full_Report_with_Intercept.xlsx\n');
catch
    writecell(metrics_summary, 'Summary_Metrics_Full_Report_with_Intercept.csv');
    fprintf('=> 成功生成严格排序排版的报表: Summary_Metrics_Full_Report_with_Intercept.csv\n');
end

% ========== 阶段三：打印各算法计算时间 ==========
fprintf('\n==============================================================\n');
fprintf('  [阶段三] 各算法计算时间统计\n');
fprintf('==============================================================\n');

if exist('history', 'var') && isfield(history, 'time')
    time_struct = history.time;
    fnames = fieldnames(time_struct);
    for i = 1:length(fnames)
        fname = fnames{i};
        if isstruct(time_struct.(fname))
            subnames = fieldnames(time_struct.(fname));
            for j = 1:length(subnames)
                fprintf('   -> %s.%s: %8.2f 秒\n', fname, subnames{j}, time_struct.(fname).(subnames{j}));
            end
        else
            fprintf('   -> %-8s: %8.2f 秒\n', upper(fname), time_struct.(fname));
        end
    end
else
    fprintf('   (提示：未在 history 中找到 time 字段)\n');
end
fprintf('--------------------------------------------------------------\n');

%% ================= 核心统计辅助函数 =================

function norm_val = normalize_positive(arr)
    min_v = min(arr); max_v = max(arr);
    if max_v == min_v, norm_val = ones(size(arr)); else, norm_val = (arr - min_v) / (max_v - min_v); end
end

function norm_val = normalize_negative(arr)
    min_v = min(arr); max_v = max(arr);
    if max_v == min_v, norm_val = ones(size(arr)); else, norm_val = (max_v - arr) / (max_v - min_v); end
end

% 【修改点】calc_metrics 增加 Intercept 输出
function [K, Intercept, R2, RMSE, NSE, MAE, Bias, KGE] = calc_metrics(sim, obs)
    valid_idx = (obs ~= -99999) & (sim ~= -99999) & ~isnan(obs) & ~isnan(sim) & ~isinf(sim);
    sim_v = sim(valid_idx); obs_v = obs(valid_idx);
    
    if length(sim_v) < 10
        K = NaN; Intercept = NaN; R2 = NaN; RMSE = NaN; NSE = NaN; MAE = NaN; Bias = NaN; KGE = NaN; return;
    end
    
    % 线性回归：观测值作 x，模拟值作 y，得到斜率 K 和截距 Intercept
    p = polyfit(obs_v, sim_v, 1);
    K = p(1);
    Intercept = p(2);
    
    coe = corrcoef(sim_v, obs_v); r = coe(2,1); R2 = r^2;
    RMSE = sqrt(mean((sim_v - obs_v).^2));
    MAE = mean(abs(sim_v - obs_v));
    Bias = mean(sim_v - obs_v); % Mean Bias (带符号)
    
    % NSE
    numerator = sum((sim_v - obs_v).^2); denominator = sum((obs_v - mean(obs_v)).^2);
    if denominator == 0, NSE = NaN; else, NSE = 1 - (numerator / denominator); end
    
    % KGE
    mu_s = mean(sim_v); mu_o = mean(obs_v);
    sigma_s = std(sim_v); sigma_o = std(obs_v);
    if mu_o == 0 || sigma_o == 0
        KGE = NaN;
    else
        alpha = sigma_s / sigma_o;
        beta = mu_s / mu_o;
        KGE = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);
    end
end

function p_matrix = make_para_matrix(p_vec)
    p_matrix = zeros(4,6); p_matrix(3,:) = p_vec; 
end

function export_data_files_standard(hh_yr_export, daily_yr, yr, alg_name, strategy_name)
    headers_hh = {'Year','Month','Day','Hour','DOY','Don','Ta','RH','VPD','Canopy height','SW','Rn','G','wind speed','CO2 concentration','GPP','measured ET(mm/h)','Rain','LAI','DataFlag','Modeled T(mm/h)','Modeled E(mm/h)','Modeled ET(mm/h)','E/ET','rsc','rss'};
    headers_daily = {'Year','Month','Day','DOY','Ta','RH','VPD','Canopy height','SW','wind speed','CO2 concentration','LAI','rsc','rss','Rn','G','GPP','measured ET(mm/d)','Rain','Modeled T(mm/d)','Modeled E(mm/d)','Modeled ET(mm/d)','E/ET'};
    
    folder_name = sprintf('Export_%d', yr);
    file_hh = fullfile(folder_name, sprintf('%d-半小时-%s-%s.xlsx', yr, alg_name, strategy_name));
    file_daily = fullfile(folder_name, sprintf('%d-天-%s-%s.xlsx', yr, alg_name, strategy_name));
    
    try
        xlswrite(file_hh, headers_hh, 'Sheet1', 'A1'); xlswrite(file_hh, hh_yr_export, 'Sheet1', 'A2');
        xlswrite(file_daily, headers_daily, 'Sheet1', 'A1'); xlswrite(file_daily, daily_yr, 'Sheet1', 'A2');
    catch
        file_hh_csv = fullfile(folder_name, sprintf('%d-半小时-%s-%s.csv', yr, alg_name, strategy_name));
        file_daily_csv = fullfile(folder_name, sprintf('%d-天-%s-%s.csv', yr, alg_name, strategy_name));
        writematrix(hh_yr_export, file_hh_csv); writematrix(daily_yr, file_daily_csv);
    end
end