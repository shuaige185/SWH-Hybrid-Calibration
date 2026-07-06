% =========================================================================
% 图2：参数敏感性、不确定性与相关性集成分析图 
% (空心透明、形色解耦终极版，新增Y5轴及固定坐标范围)
% =========================================================================
clc; close all force;

% 1. 严格加载真实数据
if ~exist('optimization_results_3504.mat', 'file')
    error('[严重错误] 找不到 optimization_results_3504.mat 文件！');
end
load('optimization_results_3504.mat');

% 原始标准物理顺序
param_names_orig = {'$b_1$', '$b_2$', '$b_3$', '$a_1$', '$g_{0}$', '$\mathrm{LightExtCoef}$'};
param_names_plain = {'b1', 'b2', 'b3', 'a1', 'g0', 'LightExtCoef'};

% （b）图专属物理顺序重排：b1, b2, g0, LightExtCoef, a1, b3
target_order_b = [1, 2, 5, 6, 4, 3]; 
param_names_b = {'$b_1$', '$b_2$', '$g_{0}$', '$\mathrm{LightExtCoef}$', '$a_1$', '$b_3$'};

% =========================================================================
% 【第一部分】：输出各单一算法的参数统计信息（用于补充材料 Table S2）
% =========================================================================
fprintf('\n========================================================\n');
fprintf('  第一部分：各单一算法最优参数统计\n');
fprintf('========================================================\n\n');

method_list = {'bayes', 'ga', 'pso', 'ssa'};
alg_display = {'BO', 'GA', 'PSO', 'SSA'};

% ------ 1.1 各单一算法最优10组参数范围（用于Table S2）------
fprintf('【Table S2 数据】各单一算法最优10组参数范围:\n');
fprintf('%-8s %-8s %-12s %-12s %-12s\n', '算法', '参数', '最小值', '最大值', '最优值');
fprintf('%-8s %-8s %-12s %-12s %-12s\n', '----', '----', '------', '------', '------');

stats_top10 = {};

for m = 1:length(method_list)
    method = method_list{m};
    disp_name = alg_display{m};
    
    all_p = history.all_methods.(method).all_params;
    all_f = history.all_methods.(method).all_fvals;
    [~, sort_idx] = sort(all_f);
    sorted_p = all_p(sort_idx, :);
    
    % 取前10组
    n_top = min(10, size(sorted_p, 1));
    top_p = sorted_p(1:n_top, :);
    best_p = top_p(1, :);
    
    for p = 1:6
        min_val = min(top_p(:, p));
        max_val = max(top_p(:, p));
        best_val = best_p(p);
        stats_top10 = [stats_top10; {disp_name, param_names_plain{p}, ...
            round(min_val, 4), round(max_val, 4), round(best_val, 4)}];
        fprintf('%-8s %-8s %12.4f %12.4f %12.4f\n', ...
            disp_name, param_names_plain{p}, min_val, max_val, best_val);
    end
    fprintf('\n');
end

% ------ 1.2 各单一算法最优10组参数统计（用于正文描述）------
fprintf('【正文描述用】各单一算法最优10组参数统计:\n');
fprintf('%-8s %-12s %-12s %-12s %-12s %-12s\n', ...
    '算法', 'b3最小值', 'b3最大值', 'b3均值', 'SGPI最小', 'SGPI最大');

for m = 1:length(method_list)
    method = method_list{m};
    disp_name = alg_display{m};
    
    all_p = history.all_methods.(method).all_params;
    all_f = history.all_methods.(method).all_fvals;
    [sorted_f, sort_idx] = sort(all_f);
    sorted_p = all_p(sort_idx, :);
    
    n_top = min(10, size(sorted_p, 1)); 
    top_p = sorted_p(1:n_top, :);
    top_f = sorted_f(1:n_top);
    
    fprintf('%-8s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
        disp_name, min(top_p(:,3)), max(top_p(:,3)), mean(top_p(:,3)), ...
        min(top_f), max(top_f));
end
fprintf('\n');

% ------ 1.3 各单一算法最优参数汇总（用于对比）------
fprintf('各单一算法最优参数汇总（目标函数值最小的一组）:\n');
fprintf('%-8s %-8s %-8s %-8s %-8s %-8s %-8s %-10s\n', ...
    '算法', 'b1', 'b2', 'b3', 'a1', 'g0', 'LightExt', 'SGPI');

for m = 1:length(method_list)
    method = method_list{m};
    disp_name = alg_display{m};
    
    all_p = history.all_methods.(method).all_params;
    all_f = history.all_methods.(method).all_fvals;
    [sorted_f, sort_idx] = sort(all_f);
    sorted_p = all_p(sort_idx, :);
    
    best_p = sorted_p(1, :);
    best_f = sorted_f(1);
    
    fprintf('%-8s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %10.6f\n', ...
        disp_name, best_p(1), best_p(2), best_p(3), best_p(4), best_p(5), best_p(6), best_f);
end
fprintf('\n');

% ------ 1.4 融合算法最优参数 ------
fprintf('【融合算法】全局最优参数:\n');
fprintf('%-8s %-8s %-8s %-8s %-8s %-8s %-8s %-10s\n', ...
    '算法', 'b1', 'b2', 'b3', 'a1', 'g0', 'LightExt', 'SGPI');
fprintf('%-8s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %10.6f\n', ...
    '融合', best_params(1), best_params(2), best_params(3), ...
    best_params(4), best_params(5), best_params(6), history.best_f);

fprintf('\n');

% =========================================================================
% 【第二部分】：输出支持正文关键数字的数据
% =========================================================================
fprintf('========================================================\n');
fprintf('  第二部分：支持正文结论的关键数字\n');
fprintf('========================================================\n\n');

fprintf('1. 各单一算法 b3 范围（基于最优10组）:\n');
for m = 1:length(method_list)
    method = method_list{m};
    disp_name = alg_display{m};
    all_p = history.all_methods.(method).all_params;
    all_f = history.all_methods.(method).all_fvals;
    [~, sort_idx] = sort(all_f);
    sorted_p = all_p(sort_idx, :);
    top_p = sorted_p(1:min(10, size(sorted_p,1)), :);
    fprintf('   %s: b3 = %.2f ~ %.2f\n', disp_name, min(top_p(:,3)), max(top_p(:,3)));
end

fprintf('\n2. 各单一算法 SGPI 范围（基于最优10组）:\n');
for m = 1:length(method_list)
    method = method_list{m};
    disp_name = alg_display{m};
    all_f = history.all_methods.(method).all_fvals;
    [sorted_f, ~] = sort(all_f);
    top_f = sorted_f(1:min(10, length(sorted_f)));
    fprintf('   %s: SGPI = %.4f ~ %.4f\n', disp_name, min(top_f), max(top_f));
end

fprintf('\n3. 融合算法各参数值:\n');
fprintf('   b1=%.2f, b2=%.2f, b3=%.2f, a1=%.2f, g0=%.4f, LightExt=%.3f, SGPI=%.4f\n', ...
    best_params(1), best_params(2), best_params(3), best_params(4), ...
    best_params(5), best_params(6), history.best_f);

% =========================================================================
% 【保存】将 Top10 统计结果保存为 Excel（Table S2 数据源）
% =========================================================================
if ~isempty(stats_top10)
    T_top10 = cell2table(stats_top10, ...
        'VariableNames', {'Algorithm', 'Parameter', 'Min', 'Max', 'Best'});
    writetable(T_top10, 'TableS2_Top10_Parameter_Ranges.xlsx');
    fprintf('\n=> Table S2 数据已保存至 TableS2_Top10_Parameter_Ranges.xlsx\n');
end
fprintf('========================================================\n\n');

% =========================================================================
% 【核心提取】：提取绘图数据 —— 4种单一算法各自前10组（总计40组）
% =========================================================================
fprintf('开始生成图2...\n');

ind_p_pool = []; 
method_list_plot = {'bayes', 'ga', 'pso', 'ssa'}; % 仅用4种单一算法绘图

for m_idx = 1:length(method_list_plot)
    m_name = method_list_plot{m_idx};
    sub_method = history.all_methods.(m_name);
    
    if isfield(sub_method, 'all_params') && ~isempty(sub_method.all_params) && ...
       isfield(sub_method, 'all_fvals') && ~isempty(sub_method.all_fvals)
        
        all_p = sub_method.all_params;
        all_f = sub_method.all_fvals;
        [~, s_idx] = sort(all_f); 
        
        n_top = min(10, length(all_f)); % ★ 严格提取每种算法前10组
        top_10 = all_p(s_idx(1:n_top), :);
        ind_p_pool = [ind_p_pool; top_10]; % 最终40组数据
    end
end

if isempty(ind_p_pool)
    error('缺失独立优化算法的真实数据！'); 
end

top_p = ind_p_pool; % 4种算法 × 10组 = 40组数据

% 动态解算真实的敏感度 SI（基于融合算法最优解）
SI = zeros(6,1);
try
    sh_cal = shuttle_carlibration; Q_val = Qs; if isempty(Q_val), Q_val = 0.46; end
    base_p = zeros(4,6); base_p(3,:) = best_params;
    [~, base_rmse] = halfhour_shuttleworth_validation(sh_cal, base_p, Q_val);
    for p_idx = 1:6
        p_pert = best_params; p_pert(p_idx) = p_pert(p_idx) * 1.1; 
        p_mat = zeros(4,6); p_mat(3,:) = p_pert;
        [~, rmse_pert] = halfhour_shuttleworth_validation(sh_cal, p_mat, Q_val);
        SI(p_idx) = abs(rmse_pert - base_rmse) / (base_rmse + 1e-12);
    end
    if sum(SI) > 0, SI = SI / sum(SI); else, error('敏感度全为0'); end
catch ME
    error('敏感度解算失败: %s', ME.message);
end

% 2. 建立高规格学术画布
fig = figure('Position', [100, 100, 1280, 960], 'Color', 'w');

% ------ 面板 (a): 参数敏感性分析 ------
ax1 = axes('Position', [0.12, 0.58, 0.33, 0.33]);
[s_SI, s_idx] = sort(SI, 'ascend');
barh(ax1, s_SI, 0.55, 'FaceColor', '#0072BD', 'EdgeColor', 'k', 'LineWidth', 1.2);
xlim(ax1, [0, max(s_SI)*1.28]); ylim(ax1, [0.5, 6.5]);
set(ax1, 'YTick', 1:6, 'YTickLabel', param_names_orig(s_idx), 'TickLabelInterpreter', 'latex', 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);
xlabel(ax1, 'Normalized Sensitivity Index (SI)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
text(ax1, 0.5, 1.08, '(a) Parameter Sensitivity Analysis', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
grid(ax1, 'on'); ax1.GridLineStyle = ':'; ax1.GridAlpha = 0.4;
for i = 1:6
    text(ax1, s_SI(i) + max(s_SI)*0.02, i, sprintf('%.3f', s_SI(i)), 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
end

% ------ 面板 (b): 5个Y轴 空心透明形色解耦排版 ------
ax2_main_pos = [0.62, 0.58, 0.24, 0.33]; 

top_p_ord = top_p(:, target_order_b); 
best_ord = best_params(target_order_b);
p_min = min(top_p_ord, [], 1); p_max = max(top_p_ord, [], 1); p_mean = mean(top_p_ord, 1);

% 五大颜色体系宣告
col_y1 = 'k';           % Y1 黑色 (b1, b2)
col_y2 = '#D95319';     % Y2 橙色 (b3)
col_y3 = '#0072BD';     % Y3 蓝色 (g0)
col_y4 = '#77AC30';     % Y4 绿色 (LightExtCoef)
col_y5 = '#7E2F8E';     % Y5 紫色 (a1)   ★ 新增

% 【第一层：Y1 主轴 (黑色，左内轴)】
ax2_y1 = axes('Position', ax2_main_pos); hold(ax2_y1, 'on');
p_y1_idx = [1, 2]; % ★ 仅保留 b1, b2，去掉了 a1
p_mean_y1 = NaN(1,6); p_min_err_y1 = NaN(1,6); p_max_err_y1 = NaN(1,6); best_p_y1 = NaN(1,6);
p_mean_y1(p_y1_idx) = p_mean(p_y1_idx); p_min_err_y1(p_y1_idx) = p_mean(p_y1_idx) - p_min(p_y1_idx); p_max_err_y1(p_y1_idx) = p_max(p_y1_idx) - p_mean(p_y1_idx); best_p_y1(p_y1_idx) = best_ord(p_y1_idx);

errorbar(ax2_y1, 1:6, p_mean_y1, p_min_err_y1, p_max_err_y1, 'o', 'Color', col_y1, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
plot(ax2_y1, 1:6, best_p_y1, 'p', 'Color', col_y1, 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
set(ax2_y1, 'YColor', col_y1, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);

% 【第二层：Y5 数据绘制层 (紫色，a1)】 ★ 新增
ax2_y5_data = axes('Position', ax2_main_pos, 'Color', 'none', 'XColor', 'none', 'YColor', 'none'); hold(ax2_y5_data, 'on');
p_y5_idx = 5; % 对应 a1
p_mean_y5 = NaN(1,6); p_min_err_y5 = NaN(1,6); p_max_err_y5 = NaN(1,6); best_p_y5 = NaN(1,6);
p_mean_y5(p_y5_idx) = p_mean(p_y5_idx); p_min_err_y5(p_y5_idx) = p_mean(p_y5_idx) - p_min(p_y5_idx); p_max_err_y5(p_y5_idx) = p_max(p_y5_idx) - p_mean(p_y5_idx); best_p_y5(p_y5_idx) = best_ord(p_y5_idx);

errorbar(ax2_y5_data, 1:6, p_mean_y5, p_min_err_y5, p_max_err_y5, 'o', 'Color', col_y5, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
plot(ax2_y5_data, 1:6, best_p_y5, 'p', 'Color', col_y5, 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

% 【第三层：Y5 视觉偏移层 (紫色，左外轴)】 ★ 新增
ax2_y5_dummy = axes('Position', [0.51, 0.58, 0.35, 0.33], 'Color', 'none', 'XColor', 'none', 'YAxisLocation', 'left', 'YColor', col_y5);
set(ax2_y5_dummy, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);

% 【第四层：Y4 数据绘制层 (绿色，LightExtCoef)】
ax2_y4_data = axes('Position', ax2_main_pos, 'Color', 'none', 'XColor', 'none', 'YColor', 'none'); hold(ax2_y4_data, 'on');
p_y4_idx = 4; % 对应 LightExtCoef
p_mean_y4 = NaN(1,6); p_min_err_y4 = NaN(1,6); p_max_err_y4 = NaN(1,6); best_p_y4 = NaN(1,6);
p_mean_y4(p_y4_idx) = p_mean(p_y4_idx); p_min_err_y4(p_y4_idx) = p_mean(p_y4_idx) - p_min(p_y4_idx); p_max_err_y4(p_y4_idx) = p_max(p_y4_idx) - p_mean(p_y4_idx); best_p_y4(p_y4_idx) = best_ord(p_y4_idx);

errorbar(ax2_y4_data, 1:6, p_mean_y4, p_min_err_y4, p_max_err_y4, 'o', 'Color', col_y4, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
plot(ax2_y4_data, 1:6, best_p_y4, 'p', 'Color', col_y4, 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

% 【第五层：Y4 视觉偏移层 (绿色，左外轴)】
ax2_y4_dummy = axes('Position', [0.57, 0.58, 0.29, 0.33], 'Color', 'none', 'XColor', 'none', 'YAxisLocation', 'left', 'YColor', col_y4);
set(ax2_y4_dummy, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);
ylabel(ax2_y4_dummy, 'Parameter Range', 'Color', 'k', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% 【第六层：Y2 右轴 (橙色，右内轴)】
ax2_y2 = axes('Position', ax2_main_pos, 'Color', 'none', 'XColor', 'none', 'YAxisLocation', 'right', 'YColor', col_y2); hold(ax2_y2, 'on');
p_y2_idx = 6; % 对应 b3
p_mean_y2 = NaN(1,6); p_min_err_y2 = NaN(1,6); p_max_err_y2 = NaN(1,6); best_p_y2 = NaN(1,6);
p_mean_y2(p_y2_idx) = p_mean(p_y2_idx); p_min_err_y2(p_y2_idx) = p_mean(p_y2_idx) - p_min(p_y2_idx); p_max_err_y2(p_y2_idx) = p_max(p_y2_idx) - p_mean(p_y2_idx); best_p_y2(p_y2_idx) = best_ord(p_y2_idx);

errorbar(ax2_y2, 1:6, p_mean_y2, p_min_err_y2, p_max_err_y2, 'o', 'Color', col_y2, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
plot(ax2_y2, 1:6, best_p_y2, 'p', 'Color', col_y2, 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
set(ax2_y2, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);

% 【第七层：Y3 数据绘制层 (蓝色，g0)】
ax2_y3_data = axes('Position', ax2_main_pos, 'Color', 'none', 'XColor', 'none', 'YColor', 'none'); hold(ax2_y3_data, 'on');
p_y3_idx = 3; % 对应 g0
p_mean_y3 = NaN(1,6); p_min_err_y3 = NaN(1,6); p_max_err_y3 = NaN(1,6); best_p_y3 = NaN(1,6);
p_mean_y3(p_y3_idx) = p_mean(p_y3_idx); p_min_err_y3(p_y3_idx) = p_mean(p_y3_idx) - p_min(p_y3_idx); p_max_err_y3(p_y3_idx) = p_max(p_y3_idx) - p_mean(p_y3_idx); best_p_y3(p_y3_idx) = best_ord(p_y3_idx);

errorbar(ax2_y3_data, 1:6, p_mean_y3, p_min_err_y3, p_max_err_y3, 'o', 'Color', col_y3, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
plot(ax2_y3_data, 1:6, best_p_y3, 'p', 'Color', col_y3, 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

% 【第八层：Y3 视觉偏移层 (蓝色，右外轴)】
ax2_y3_dummy = axes('Position', [0.62, 0.58, 0.29, 0.33], 'Color', 'none', 'XColor', 'none', 'YAxisLocation', 'right', 'YColor', col_y3);
set(ax2_y3_dummy, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);

% --- 同步所有轴的 X 与 Y limits ---
xlim(ax2_y1, [0.5, 6.5]); xlim(ax2_y2, [0.5, 6.5]); 
xlim(ax2_y3_data, [0.5, 6.5]); xlim(ax2_y3_dummy, [0.5, 6.5]);
xlim(ax2_y4_data, [0.5, 6.5]); xlim(ax2_y4_dummy, [0.5, 6.5]);
xlim(ax2_y5_data, [0.5, 6.5]); xlim(ax2_y5_dummy, [0.5, 6.5]);

% 【严格按照用户要求设定每个Y轴的起点和终点】
% Y1 (b1, b2) 范围: 1.00-5.00
ylim(ax2_y1, [1.0, 5.0]);

% Y2 (b3) 范围: 1.00-1000.00
ylim(ax2_y2, [1.0, 1000.0]);

% Y3 (g0) 范围: 0.001-0.100
ylim(ax2_y3_data, [0.001, 0.100]); ylim(ax2_y3_dummy, [0.001, 0.100]);

% Y4 (LightExtCoef) 范围: 0.450-0.700
ylim(ax2_y4_data, [0.450, 0.700]); ylim(ax2_y4_dummy, [0.450, 0.700]);

% Y5 (a1) 范围: 1.00-100.00
ylim(ax2_y5_data, [1.0, 100.0]); ylim(ax2_y5_dummy, [1.0, 100.0]);

% 彻底禁用自带的 X 轴标签
set(ax2_y1, 'XTick', 1:6, 'XTickLabel', []); 
grid(ax2_y1, 'on'); ax2_y1.GridLineStyle = ':'; ax2_y1.GridAlpha = 0.5;

% 【彩色 X 轴标签黑科技，a1标签变为Y5的颜色】
y_lim = ylim(ax2_y1);
y_pos = y_lim(1) - (y_lim(2)-y_lim(1))*0.03; 
% 对应顺序: b1(Y1黑), b2(Y1黑), g0(Y3蓝), LightExt(Y4绿), a1(Y5紫), b3(Y2橙)
colors_b = {col_y1, col_y1, col_y3, col_y4, col_y5, col_y2}; 

for i = 1:6
    text(ax2_y1, i, y_pos, param_names_b{i}, 'Color', colors_b{i}, ...
        'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman', ...
        'Rotation', 45, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Clipping', 'off');
end

text(ax2_y1, 0.5, 1.08, '(b) Parameter Uncertainty Analysis', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% 【中立空心图注】
h_err = errorbar(ax2_y1, nan, nan, nan, 'o', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
h_opt = plot(ax2_y1, nan, nan, 'p', 'Color', [0.4 0.4 0.4], 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
lgd = legend(ax2_y1, [h_err, h_opt], {'Ensemble Range', 'Global Optimal'}, ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', 13, 'FontName', 'Times New Roman');
drawnow; 
lgd.Position(1) = lgd.Position(1) + 0.10; 

% ------ 面板 (c): 参数相关性分析 ------
ax3 = axes('Position', [0.12, 0.11, 0.33, 0.33]);
R = corrcoef(top_p(:, 1:6));
for k=1:6, R(k,k) = 1.00; end 
imagesc(ax3, R); colormap(ax3, 'parula'); 
cb = colorbar(ax3); set(cb, 'FontSize', 13, 'FontName', 'Times New Roman');
set(ax3, 'XTick', 1:6, 'XTickLabel', param_names_orig, 'YTick', 1:6, 'YTickLabel', param_names_orig);
set(ax3, 'TickLabelInterpreter', 'latex', 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);
xtickangle(ax3, 45); 
text(ax3, 0.5, 1.08, '(c) Parameter Correlation Analysis', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

for i = 1:6
    for j = 1:6
        txt_col = 'k'; if abs(R(i,j)) > 0.55, txt_col = 'w'; end
        if i == j
            text(ax3, j, i, '1.00', 'HorizontalAlignment', 'center', 'Color', 'k', 'BackgroundColor', [1 1 1 0.6], 'FontWeight', 'bold', 'FontSize', 12);
        else
            text(ax3, j, i, sprintf('%.2f', R(i,j)), 'HorizontalAlignment', 'center', 'Color', txt_col, 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
end

% =========================================================================
% 导出图表（自动保存至当前文件夹）
% =========================================================================
save_dir = pwd; 
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
save_path = fullfile(save_dir, 'Figure_2_Integrated_Uncertainty.emf');
print(fig, save_path, '-dmeta', '-r600');
fprintf('=> 图2 已成功保存至: %s\n', save_path);

fprintf('\n========================================================\n');
fprintf('  全部完成！\n');
fprintf('  输出文件：\n');
fprintf('  1. Figure_2_Integrated_Uncertainty.emf (图2)\n');
fprintf('  2. TableS2_Top10_Parameter_Ranges.xlsx (补充材料Table S2)\n');
fprintf('========================================================\n');