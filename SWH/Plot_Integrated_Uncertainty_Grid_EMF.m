% =========================================================================
% 图2：参数敏感性、不确定性与相关性集成分析图 (空心透明、形色解耦终极版)
% =========================================================================
clc; close all force;

% 1. 严格加载真实数据
if ~exist('optimization_results_3504.mat', 'file')
    error('[严重错误] 找不到 optimization_results_3504.mat 文件！');
end
load('optimization_results_3504.mat');

% 原始标准物理顺序
param_names_orig = {'$b_1$', '$b_2$', '$b_3$', '$a_1$', '$g_{0}$', '$\mathrm{LightExtCoef}$'};

% （b）图专属物理顺序重排：b1, b2, g0, LightExtCoef, a1, b3
target_order_b = [1, 2, 5, 6, 4, 3]; 
param_names_b = {'$b_1$', '$b_2$', '$g_{0}$', '$\mathrm{LightExtCoef}$', '$a_1$', '$b_3$'};

% 提取全局大池子高质量解
all_p_pool = [];    
all_fvals_pool = []; 

if exist('history', 'var') && isfield(history, 'all_methods')
    method_names = fieldnames(history.all_methods);
    for m_idx = 1:length(method_names)
        m_name = method_names{m_idx};
        sub_method = history.all_methods.(m_name);
        if isfield(sub_method, 'all_params') && ~isempty(sub_method.all_params) && ...
           isfield(sub_method, 'all_fvals') && ~isempty(sub_method.all_fvals)
            all_p_pool = [all_p_pool; sub_method.all_params];
            all_fvals_pool = [all_fvals_pool; sub_method.all_fvals];
        end
    end
end
if isempty(all_p_pool), error('缺失真实优化数据！'); end

% 严格筛选表现最好的前 100 名
[~, sort_idx] = sort(all_fvals_pool, 'ascend');
top_n = min(100, length(sort_idx));
top_p = all_p_pool(sort_idx(1:top_n), :); 

% 动态解算真实的敏感度 SI 
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


% ------ 面板 (b): 4个Y轴 空心透明形色解耦排版 ------
% 主绘图区
ax2_main_pos = [0.62, 0.58, 0.24, 0.33]; 

top_p_ord = top_p(:, target_order_b); 
best_ord = best_params(target_order_b);
p_min = min(top_p_ord, [], 1); p_max = max(top_p_ord, [], 1); p_mean = mean(top_p_ord, 1);

% 四大颜色体系宣告
col_y1 = 'k';           % Y1 黑色 (b1, b2, a1)
col_y2 = '#D95319';     % Y2 橙色 (b3)
col_y3 = '#0072BD';     % Y3 蓝色 (g0)
col_y4 = '#77AC30';     % Y4 绿色 (LightExtCoef)

% 【第一层：Y1 主轴 (黑色，左内轴)】
ax2_y1 = axes('Position', ax2_main_pos); hold(ax2_y1, 'on');
p_y1_idx = [1, 2, 5]; % 对应 b1, b2, a1
p_mean_y1 = NaN(1,6); p_min_err_y1 = NaN(1,6); p_max_err_y1 = NaN(1,6); best_p_y1 = NaN(1,6);
p_mean_y1(p_y1_idx) = p_mean(p_y1_idx); p_min_err_y1(p_y1_idx) = p_mean(p_y1_idx) - p_min(p_y1_idx); p_max_err_y1(p_y1_idx) = p_max(p_y1_idx) - p_mean(p_y1_idx); best_p_y1(p_y1_idx) = best_ord(p_y1_idx);

% 空心化：MarkerFaceColor 设为 'none'，LineWidth 加粗为 1.5 增加质感
errorbar(ax2_y1, 1:6, p_mean_y1, p_min_err_y1, p_max_err_y1, 'o', 'Color', col_y1, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
plot(ax2_y1, 1:6, best_p_y1, 'p', 'Color', col_y1, 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
set(ax2_y1, 'YColor', col_y1, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);

% 【第二层：Y4 数据绘制层 (绿色，LightExtCoef)】
ax2_y4_data = axes('Position', ax2_main_pos, 'Color', 'none', 'XColor', 'none', 'YColor', 'none'); hold(ax2_y4_data, 'on');
p_y4_idx = 4; % 对应 LightExtCoef
p_mean_y4 = NaN(1,6); p_min_err_y4 = NaN(1,6); p_max_err_y4 = NaN(1,6); best_p_y4 = NaN(1,6);
p_mean_y4(p_y4_idx) = p_mean(p_y4_idx); p_min_err_y4(p_y4_idx) = p_mean(p_y4_idx) - p_min(p_y4_idx); p_max_err_y4(p_y4_idx) = p_max(p_y4_idx) - p_mean(p_y4_idx); best_p_y4(p_y4_idx) = best_ord(p_y4_idx);

errorbar(ax2_y4_data, 1:6, p_mean_y4, p_min_err_y4, p_max_err_y4, 'o', 'Color', col_y4, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
plot(ax2_y4_data, 1:6, best_p_y4, 'p', 'Color', col_y4, 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

% 【第三层：Y4 视觉偏移层 (绿色，左外轴)】
ax2_y4_dummy = axes('Position', [0.57, 0.58, 0.29, 0.33], 'Color', 'none', 'XColor', 'none', 'YAxisLocation', 'left', 'YColor', col_y4);
set(ax2_y4_dummy, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);
ylabel(ax2_y4_dummy, 'Parameter Range', 'Color', 'k', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% 【第四层：Y2 右轴 (橙色，右内轴)】
ax2_y2 = axes('Position', ax2_main_pos, 'Color', 'none', 'XColor', 'none', 'YAxisLocation', 'right', 'YColor', col_y2); hold(ax2_y2, 'on');
p_y2_idx = 6; % 对应 b3
p_mean_y2 = NaN(1,6); p_min_err_y2 = NaN(1,6); p_max_err_y2 = NaN(1,6); best_p_y2 = NaN(1,6);
p_mean_y2(p_y2_idx) = p_mean(p_y2_idx); p_min_err_y2(p_y2_idx) = p_mean(p_y2_idx) - p_min(p_y2_idx); p_max_err_y2(p_y2_idx) = p_max(p_y2_idx) - p_mean(p_y2_idx); best_p_y2(p_y2_idx) = best_ord(p_y2_idx);

errorbar(ax2_y2, 1:6, p_mean_y2, p_min_err_y2, p_max_err_y2, 'o', 'Color', col_y2, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
plot(ax2_y2, 1:6, best_p_y2, 'p', 'Color', col_y2, 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
set(ax2_y2, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);

% 【第五层：Y3 数据绘制层 (蓝色，g0)】
ax2_y3_data = axes('Position', ax2_main_pos, 'Color', 'none', 'XColor', 'none', 'YColor', 'none'); hold(ax2_y3_data, 'on');
p_y3_idx = 3; % 对应 g0
p_mean_y3 = NaN(1,6); p_min_err_y3 = NaN(1,6); p_max_err_y3 = NaN(1,6); best_p_y3 = NaN(1,6);
p_mean_y3(p_y3_idx) = p_mean(p_y3_idx); p_min_err_y3(p_y3_idx) = p_mean(p_y3_idx) - p_min(p_y3_idx); p_max_err_y3(p_y3_idx) = p_max(p_y3_idx) - p_mean(p_y3_idx); best_p_y3(p_y3_idx) = best_ord(p_y3_idx);

errorbar(ax2_y3_data, 1:6, p_mean_y3, p_min_err_y3, p_max_err_y3, 'o', 'Color', col_y3, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
plot(ax2_y3_data, 1:6, best_p_y3, 'p', 'Color', col_y3, 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

% 【第六层：Y3 视觉偏移层 (蓝色，右外轴)】
ax2_y3_dummy = axes('Position', [0.62, 0.58, 0.29, 0.33], 'Color', 'none', 'XColor', 'none', 'YAxisLocation', 'right', 'YColor', col_y3);
set(ax2_y3_dummy, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);


% --- 同步所有轴的 X 与 Y limits ---
xlim(ax2_y1, [0.5, 6.5]); xlim(ax2_y2, [0.5, 6.5]); 
xlim(ax2_y3_data, [0.5, 6.5]); xlim(ax2_y3_dummy, [0.5, 6.5]);
xlim(ax2_y4_data, [0.5, 6.5]); xlim(ax2_y4_dummy, [0.5, 6.5]);

% 同步 Y3 (蓝色) 轴上下限
y3_lims = [min([p_min(3), best_ord(3)]) * 0.9, max([p_max(3), best_ord(3)]) * 1.1];
if y3_lims(2) <= y3_lims(1), y3_lims = [0, 0.1]; end
ylim(ax2_y3_data, y3_lims); ylim(ax2_y3_dummy, y3_lims);

% 同步 Y4 (绿色) 轴上下限
y4_lims = [min([p_min(4), best_ord(4)]) * 0.95, max([p_max(4), best_ord(4)]) * 1.05];
if y4_lims(2) <= y4_lims(1), y4_lims = [0.4, 0.8]; end
ylim(ax2_y4_data, y4_lims); ylim(ax2_y4_dummy, y4_lims);

% 彻底禁用自带的 X 轴标签
set(ax2_y1, 'XTick', 1:6, 'XTickLabel', []); 
grid(ax2_y1, 'on'); ax2_y1.GridLineStyle = ':'; ax2_y1.GridAlpha = 0.5;

% 【彩色 X 轴标签黑科技】
y_lim = ylim(ax2_y1);
y_pos = y_lim(1) - (y_lim(2)-y_lim(1))*0.03; 
% 对应顺序: b1(黑), b2(黑), g0(蓝), Light(绿), a1(黑), b3(橙)
colors_b = {col_y1, col_y1, col_y3, col_y4, col_y1, col_y2}; 

for i = 1:6
    text(ax2_y1, i, y_pos, param_names_b{i}, 'Color', colors_b{i}, ...
        'Interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times New Roman', ...
        'Rotation', 45, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Clipping', 'off');
end

text(ax2_y1, 0.5, 1.08, '(b) Parameter Uncertainty Analysis', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% 【中立空心图注，并向右平移防遮挡】
% 使用深灰色 [0.4 0.4 0.4] 彻底消除偏袒某一种颜色的误导
h_err = errorbar(ax2_y1, nan, nan, nan, 'o', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'none');
h_opt = plot(ax2_y1, nan, nan, 'p', 'Color', [0.4 0.4 0.4], 'MarkerSize', 11, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

lgd = legend(ax2_y1, [h_err, h_opt], {'Ensemble Range', 'Global Optimal'}, ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', 13, 'FontName', 'Times New Roman');
drawnow; 
% 强行将图注框向右平移，释放左侧空间，完全不遮挡！
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

% 设置保存目录（修改为当前目录下的 Figures 文件夹）
save_dir = fullfile(pwd, 'Figures');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% 拼接完整路径（注意这里要改成图2的名字，而不是图1）
save_path = fullfile(save_dir, 'Figure_2_Integrated_Uncertainty.emf');

% 使用 -dmeta 参数导出 EMF（MATLAB 2019b 完美兼容，最适合插入 Word）
print(fig, save_path, '-dmeta', '-r600');
fprintf('=> [图2] 已成功保存至: %s\n', save_path);