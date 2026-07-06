% =========================================================================
% 图1：真实收敛曲线与运行时间（精细微调版——统一15号字、标签居中、对称轴名）
% =========================================================================
clc; close all;
try cd(fileparts(mfilename('fullpath'))); catch; end

% 加载优化核心数据
load('optimization_results_3504.mat', 'history');

% 【布局微调】：由于字号改为15号，将画布缩减至最适合论文排版的紧凑尺寸，防止画面太空旷
fig = figure('Position', [100, 100, 1000, 480], 'Color', 'w');

% ------ 面板 (a): 收敛曲线 (绝对坐标控制) ------
ax1 = axes('Position', [0.11, 0.16, 0.37, 0.71]); 
hold on;

colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'}; 
line_styles = {'-', '--', '-.', ':'}; 
disp_names = {'BO', 'GA', 'PSO', 'SSA'};
methods = {'bayes', 'ga', 'pso', 'ssa'};
max_y = 0; 

for i = 1:4
    m = methods{i};
    if isfield(history, 'all_methods') && isfield(history.all_methods, m)
        fvals = history.all_methods.(m).all_fvals;
        best_fvals = zeros(length(fvals), 1); 
        best_fvals(1) = fvals(1);
        for j = 2:length(fvals)
            best_fvals(j) = min(best_fvals(j-1), fvals(j)); 
        end
        plot(1:length(fvals), best_fvals, 'Color', colors{i}, 'LineStyle', line_styles{i}, ...
            'LineWidth', 2.5, 'DisplayName', disp_names{i});
        max_y = max(max_y, max(best_fvals));
    end
end
hold off;

% 严格对数坐标系与范围
set(ax1, 'XScale', 'log'); 
xlim(ax1, [1, 50000]);

% 规范化对数刻度
log_ticks = [1, 10, 100, 1000, 10000, 50000];
set(ax1, 'XTick', log_ticks);
set(ax1, 'XTickLabel', {'10^0', '10^1', '10^2', '10^3', '10^4', '5\times10^4'});

% 纵坐标刻度
y_max_limit = ceil(max_y * 1.15);
if y_max_limit < 3, y_max_limit = 3; end 
set(ax1, 'YTick', 0:1:y_max_limit);
set(ax1, 'YMinorTick', 'off'); 
ylim(ax1, [0, y_max_limit]);

% 【字体修改】：坐标轴数字字号统一改为 15 号新罗马体
set(ax1, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);

% 【字体修改】：标签字号统一改为 15 号加粗
xlabel(ax1, 'Evaluations', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel(ax1, '$\mathbf{SGPI_{obj}}$', 'Interpreter', 'latex', 'FontSize', 15);
% 【重点修改】：保持原垂直高度(1.08)，将水平坐标设为0.5并设置居中对齐，字号统一改为 15 号
text(ax1, 0.5, 1.08, '(a)', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% 【字体修改】：图例字号统一改为 15 号
legend(ax1, 'Location', 'northeast', 'Box', 'off', 'FontSize', 15, 'FontName', 'Times New Roman');
grid(ax1, 'on'); ax1.GridLineStyle = ':'; ax1.GridAlpha = 0.5;


% ------ 面板 (b): 运行时间条形图 (绝对坐标控制) ------
ax2 = axes('Position', [0.59, 0.16, 0.37, 0.71]); 

time_data = zeros(4, 1);
if isfield(history, 'time')
    if isfield(history.time, 'bayes'), time_data(1) = history.time.bayes; end
    if isfield(history.time, 'ga'), time_data(2) = history.time.ga; end
    if isfield(history.time, 'pso'), time_data(3) = history.time.pso; end
    if isfield(history.time, 'ssa'), time_data(4) = history.time.ssa; end
end

% 绘制条形图
b = bar(ax2, 1:4, time_data, 0.45, 'FaceColor', '#0072BD', 'EdgeColor', 'k', 'LineWidth', 1.2);

set(ax2, 'XTick', 1:4, 'XTickLabel', disp_names);
ylim(ax2, [0, max(max(time_data)*1.25, 10)]);

% 【字体修改】：右图坐标轴数字字号统一改为 15 号
set(ax2, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);

% 【重点修改】：在(b)图下方正式添加英文 X 轴标签 'Algorithms'，与左图对称，字号 15 号加粗
xlabel(ax2, 'Algorithms', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel(ax2, 'Execution Time (s)', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman'); 

% 【重点修改】：保持原垂直高度(1.08)，将水平坐标设为0.5并设置居中对齐，字号统一改为 15 号
text(ax2, 0.5, 1.08, '(b)', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

grid(ax2, 'on'); ax2.GridLineStyle = ':'; ax2.GridAlpha = 0.3;

% 【字体修改】：柱状图顶部的具体数值字号统一改为 15 号
for i = 1:4
    if time_data(i) > 0
        text(ax2, i, time_data(i) + max(time_data)*0.02, sprintf('%.1f', time_data(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold');
    end
end

% =========================================================================
% 导出图表 (保存到 F:\pic)
% =========================================================================
save_dir = 'F:\pic';

% 检查文件夹是否存在，如果不存在则自动创建
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% 拼接完整路径
save_path = fullfile(save_dir, 'Figure_1_Convergence_and_Time.emf');

% 使用 -dmeta 参数导出 EMF（MATLAB 2019b 完美兼容，最适合插入 Word）
print(fig, save_path, '-dmeta', '-r600');
fprintf('=> 完美微调版 [图1] 已成功保存至: %s\n', save_path);