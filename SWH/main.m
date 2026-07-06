% Inputdata（20 colums,
% 1 year    2 month    3 day    4 hour    5 doy (day of year)    6 don (day or night, day is 1, night is -99999)
% 7 Ta (air temperature, ℃)    8 RH(relative humidity, %)      9 VPD(kpa)
% 10 canopy height(m)           11 SW (soil volumum water content m3m-3)
% 12 Rn (net raidation, Wa m-2)     13 G (soil heat flux, Wa m-2)    14 windspeed (m s-1) 
% 15 co2( mg Co2 m-3),can set it as 500          16 GPP(mg CO2 m-2s-1)         17 ET_measure (gH2O m-2s-1)
% 18 Rain (mm),not mandatory                 19 LAI             20 DataFlag (set it as 0)
% main.m 

clear all; close all; clc;

fprintf('=========================================\n');
fprintf('  约束增强多算法融合全局优化系统\n');
fprintf('  训练数据优化\n');
fprintf('=========================================\n\n');

% 用户输入（现已强制自动化）
% mode = input('Is the time step of your inputdata hourly (1) or daily (2)?\n ');
mode = 1; % 【修改点】：默认强制按小时步长处理，如需日步长可手动改为2

% ========== 读取数据 ==========
if ~exist('inputdata.txt', 'file')
    error('【严重错误】在当前目录下找不到 inputdata.txt 输入数据文件！');
end
[shuttle_original] = textread('inputdata.txt', '', 'headerlines', 1);
disp(['原始数据行数: ' num2str(size(shuttle_original, 1))]);
disp(['原始数据列数: ' num2str(size(shuttle_original, 2))]);

% 数据统计信息
fprintf('\n数据统计信息:\n');
fprintf('  总数据量: %d\n', size(shuttle_original, 1));
fprintf('  变量数量: %d\n', size(shuttle_original, 2));

if size(shuttle_original, 1) > 1
    year_start = shuttle_original(1,1);
    year_end = shuttle_original(end,1);
    doy_start = shuttle_original(1,5);
    doy_end = shuttle_original(end,5);
    fprintf('  数据时间跨度: 第%d年第%d天 到 第%d年第%d天\n', ...
        year_start, doy_start, year_end, doy_end);
end

Qs = 0.46;

% 异常数据剔除
fprintf('\n1. 数据筛选...\n');
[shuttle_backup, shuttle_FlagBadData] = filtering(shuttle_original);
fprintf('  筛选后数据量: %d (剔除了 %d 个异常点)\n', ...
    size(shuttle_backup, 1), size(shuttle_original, 1) - size(shuttle_backup, 1));

% ========== 数据分割 ==========
n_total = size(shuttle_backup, 1);
n_train = floor(n_total / 2);
shuttle_carlibration = shuttle_backup(1:n_train, :);
shuttle_validation = shuttle_backup(n_train+1:end, :);
fprintf('  训练集: %d 组数据\n', size(shuttle_carlibration, 1));
fprintf('  验证集: %d 组数据\n', size(shuttle_validation, 1));

% ========== 参数优化 ==========
fprintf('\n3. 开始参数优化...\n');
fprintf('   使用约束增强多算法融合优化\n');

% 设置优化超参数
max_iterations = 250;  % 总迭代次数
constraint_weight = 0;  % 约束权重

% 开始计时
tic;

% 调用优化函数
[best_params, para_range_10, history, opt_info, all_method_params] = ...
    bayesian_optimization_with_global_search_constrained_3504(...
        shuttle_carlibration, Qs, max_iterations, constraint_weight);
% 结束计时
elapsed_time = toc;
fprintf('\n优化完成，耗时 %.2f 分钟\n', elapsed_time/60);

% 显示各阶段时间
fprintf('\n========== 各算法运行时间 ==========\n');
if isfield(history, 'time')
    fnames = fieldnames(history.time);
    for i = 1:length(fnames)
        if isstruct(history.time.(fnames{i}))
            subnames = fieldnames(history.time.(fnames{i}));
            for j = 1:length(subnames)
                fprintf('  %s.%s: %.2f 秒\n', fnames{i}, subnames{j}, history.time.(fnames{i}).(subnames{j}));
            end
        else
            fprintf('  %s: %.2f 秒\n', fnames{i}, history.time.(fnames{i}));
        end
    end
    fprintf('  融合优化总耗时: %.2f 秒 (%.2f 分钟)\n', elapsed_time, elapsed_time/60);
else
    fprintf('  未记录时间信息\n');
end

% 检查优化结果
if ~exist('best_params', 'var') || isempty(best_params) || ...
   ~exist('para_range_10', 'var') || isempty(para_range_10)
    fprintf('\n警告：优化返回空值，但已从输出中获取参数\n');
    if exist('history', 'var') && isfield(history, 'best_x')
        best_params = history.best_x;
        fprintf('从history中获取best_params\n');
    end
    if ~exist('para_range_10', 'var') || isempty(para_range_10)
        para_range_10 = zeros(4, 6);
        para_range_10(3, :) = best_params;
        para_range_10(1, :) = best_params * 0.5;
        para_range_10(2, :) = best_params * 2.0;
        para_range_10(4, :) = best_params * 0.2;
    end
end

% ========== 显示可行性分析结果 ==========
fprintf('\n========== 可行性分析结果 ==========\n');
if isfield(history, 'feasibility_summary')
    method_names = fieldnames(history.feasibility_summary);
    fprintf('各方法fmincon优化参数的可行性:\n');
    for i = 1:length(method_names)
        method_name = method_names{i};
        summary = history.feasibility_summary.(method_name);
        if summary.is_feasible
            fprintf('  %s: ? 可行 (违反0个约束)\n', method_name);
        else
            fprintf('  %s: ? 不可行 (违反%d个约束, 违反量: %.4f)\n', ...
                method_name, summary.constraint_violations, summary.violation_magnitude);
        end
    end
end

if isfield(history, 'is_valid')
    if history.is_valid
        fprintf('\n? 最佳参数满足所有约束条件\n');
    else
        fprintf('\n? 注意：最佳参数不完全满足所有约束条件\n');
    end
end

% ========== 模型验证和计算 ==========
fprintf('\n4. 开始模型验证和计算...\n');
[shuttle_validation_result, RMSE] = halfhour_shuttleworth_validation(shuttle_validation, para_range_10, Qs);
[shuttle_halfhour] = halfhour_shuttleworth_fillgap(shuttle_FlagBadData, para_range_10, Qs);

% ========== 保存结果 ==========
fprintf('\n5. 保存结果...\n');
try
    save('optimization_results_3504.mat', 'best_params', 'para_range_10', ...
        'history', 'opt_info', 'RMSE', 'shuttle_carlibration', 'shuttle_validation', ...
        'shuttle_FlagBadData', 'shuttle_original', 'Qs', 'all_method_params', '-mat');
    fprintf('   ? 优化结果及所有画图依赖数据已保存到 optimization_results_3504.mat\n');
    
    save('model_results_3504.mat', 'shuttle_validation_result', ...
        'shuttle_halfhour', '-mat');
    fprintf('   ? 模型结果已保存到 model_results_3504.mat\n');
    
    dlmwrite('best_params_3504.txt', best_params, 'delimiter', '\t');
    fprintf('   ? 最佳参数已保存到 best_params_3504.txt\n');
catch ME
    fprintf('   ? 警告：保存文件时出错: %s\n', ME.message);
end

% ========== 输出结果摘要 ==========
fprintf('\n========== 优化结果摘要 ==========\n');
fprintf('训练数据量: %d\n', size(shuttle_carlibration, 1));
fprintf('验证数据量: %d\n', size(shuttle_validation, 1));
fprintf('最佳目标函数值: %.4f\n', history.best_f);
if exist('opt_info', 'var') && isfield(opt_info, 'best_method')
    fprintf('最佳参数来源: %s\n', opt_info.best_method);
end

fprintf('\n最佳参数值:\n');
param_names = {'b1', 'b2', 'b3', 'a1', 'g01', 'LightExtCoef'};
for i = 1:length(best_params)
    fprintf('  %-15s = %.4f (范围: %.4f - %.4f)\n', ...
        param_names{i}, best_params(i), para_range_10(1,i), para_range_10(2,i));
end

if isfield(history, 'final_validation')
    fprintf('\n最终参数验证结果:\n');
    val = history.final_validation;
    fprintf('  ET斜率: %.3f (要求0.9-1.1)\n', val.slope_ET);
    fprintf('  ET-R?: %.3f (要求>=0.8)\n', val.R2_ET);
    fprintf('  ET-RMSE: %.4f\n', val.RMSE_ET);
    fprintf('  ET-NSE: %.3f (要求>=0.8)\n', val.NSE_ET);
    fprintf('  有效样本数: %d\n', val.ET_samples);
end
fprintf('\n===================================\n');

% ========== 6. 显示8组参数统计指标 ==========
if exist('all_method_params', 'var')
    display_8parameter_stats(all_method_params, shuttle_carlibration, Qs);
end

% ========== 7. 运行8个参数集的模型并打包 ==========
fprintf('\n7. 运行8个参数集的模型...\n');
if exist('all_method_params', 'var')
    all_simulation_results = struct(); 
    method_list = {'bayes', 'ga', 'pso', 'ssa'};
    param_suffixes = {'', '_f'};  
    success_count = 0;
    
    for m_idx = 1:length(method_list)
        for s_idx = 1:length(param_suffixes)
            param_name = [method_list{m_idx} param_suffixes{s_idx}];
            if isfield(all_method_params, param_name) && ~isempty(all_method_params.(param_name))
                fprintf('\n运行模型: %s...\n', param_name);
                try
                    params = all_method_params.(param_name);
                    para_range_temp = zeros(4, 6);
                    para_range_temp(1:3, :) = repmat(params, 3, 1);
                    
                    [shuttle_halfhour_temp] = halfhour_shuttleworth_fillgap(shuttle_FlagBadData, para_range_temp, Qs);
                    
                    try
                        [shuttledaily_temp, ~] = shuttlesum(shuttle_halfhour_temp);
                    catch ME_inner
                        fprintf('  shuttlesum运行失败: %s\n', ME_inner.message);
                        shuttledaily_temp = [];
                    end
                    
                    if ~isempty(shuttle_halfhour_temp)
                        all_simulation_results.(param_name).halfhour = shuttle_halfhour_temp;
                    end
                    if ~isempty(shuttledaily_temp)
                        all_simulation_results.(param_name).daily = shuttledaily_temp;
                    end
                    fprintf('  %s 结果已成功打包存入内存结构体\n', param_name);
                    success_count = success_count + 1;                
                catch ME
                    fprintf('  运行失败: %s\n', ME.message);
                end
            end
        end
    end
    fprintf('\n完成 %d/8 个参数集的运行\n', success_count);
    
    fprintf('正在将所有尺度模拟结果打包为一个大文件...\n');
    save('All_Simulation_Results_3504.mat', 'all_simulation_results', '-v7.3');
    fprintf(' => 成功！所有结果已合并保存为 All_Simulation_Results_3504.mat\n');
end

% ========== 8. 生成 Rebuttal 多维度终极报表 (Bias, KGE, GPI) ==========
fprintf('\n======================================================\n');
fprintf('8. 开始生成多尺度终极统计报表...\n');
try
    generate_rebuttal_exports_and_metrics;
    fprintf(' -> 恭喜！多尺度终极报表导出已成功完成！\n');
catch ME
    fprintf('   -> 生成报表时出错: %s\n', ME.message);
end

% ========== 9. 集中生成 Rebuttal 专用高清 EMF 可视化结果 ==========
fprintf('\n======================================================\n');
fprintf('9. 集中生成 Rebuttal 专用高清 EMF 矢量图...\n');
try
    if exist('Plot_Convergence_Time_EMF.m', 'file'), Plot_Convergence_Time_EMF; end
    if exist('Plot_Integrated_Uncertainty_Grid_EMF.m', 'file'), Plot_Integrated_Uncertainty_Grid_EMF; end
    if exist('Plot_SGPI_Weight_Sensitivity_EMF.m', 'file'), Plot_SGPI_Weight_Sensitivity_EMF; end
    fprintf('\n -> 完美！3张 Rebuttal 核心防守矢量图表已全部生成！\n');
catch ME
    fprintf('   -> 可视化过程中出错: %s\n', ME.message);
end

fprintf('\n======================================================\n');
fprintf('  ? 全部工作流执行完毕！请检查当前文件夹中的生成的图片和Excel表格！\n');
fprintf('======================================================\n');