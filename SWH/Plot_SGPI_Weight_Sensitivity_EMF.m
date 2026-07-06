% =========================================================================
% 图3：动态数据 SGPI 权重验证图 (绝对坐标、强力锁死刻度、专属颜色版)
% =========================================================================
clc; close all force; % 强制关闭所有旧窗口，防止干扰

% 检查数据环境
if ~exist('optimization_results_3504.mat', 'file')
    error('[严重错误] 找不到 optimization_results_3504.mat 文件！');
end
load('optimization_results_3504.mat');

if ~exist('all_method_params', 'var') || ~exist('shuttle_carlibration', 'var')
    error('[严重错误] 无法加载 all_method_params 或 shuttle_carlibration！');
end
if ~exist('Qs', 'var'), Qs = 0.46; end

fprintf('=> 正在重新计算 SGPI 动态权重数据，请稍候...\n');

scenarios = {
    '1:1:1:1', [1.0, 1.0, 1.0, 1.0];
    '1:1:3:3', [1.0, 1.0, 3.0, 3.0];
    '3:3:1:1', [3.0, 3.0, 1.0, 1.0];
    '1:3:3:1', [1.0, 3.0, 3.0, 1.0];
    '3:1:3:1', [3.0, 1.0, 3.0, 1.0]
};
num_scenarios = size(scenarios, 1);
scenario_labels = scenarios(:, 1);

method_list = {'bayes', 'ga', 'pso', 'ssa'};
valid_methods = [method_list, cellfun(@(x) [x '_f'], method_list, 'UniformOutput', false)];
num_methods = length(valid_methods);
raw_metrics = zeros(num_methods, 4);

% 提取结果数据
for m = 1:num_methods
    p_name = valid_methods{m};
    if isfield(all_method_params, p_name) && ~isempty(all_method_params.(p_name))
        p_matrix = zeros(4,6); p_matrix(3,:) = all_method_params.(p_name);
        [shuttle_out, ~] = halfhour_shuttleworth_validation(shuttle_carlibration, p_matrix, Qs);
        ET_obs = shuttle_carlibration(:, 17);
        ET_sim = shuttle_out(:, 23);
        
        valid_idx = (ET_obs ~= -99999) & (ET_sim ~= -99999) & ~isnan(ET_obs) & ~isnan(ET_sim);
        if sum(valid_idx) >= 50
            obs = ET_obs(valid_idx); sim = ET_sim(valid_idx);
            cc = corrcoef(sim, obs); if numel(cc) >= 4, raw_metrics(m, 1) = cc(2,1)^2; end
            denom = sum((obs - mean(obs)).^2);
            if denom ~= 0, raw_metrics(m, 2) = 1 - sum((sim - obs).^2) / denom; end
            raw_metrics(m, 3) = sqrt(mean((sim - obs).^2));
            raw_metrics(m, 4) = mean(abs(sim - obs));
        end
    end
end

% 归一化与打分排名
scaled_metrics = zeros(num_methods, 4);
for j = 1:4
    med = median(raw_metrics(:, j));
    mad_val = median(abs(raw_metrics(:, j) - med));
    if mad_val == 0, mad_val = std(raw_metrics(:, j)) + 1e-6; end
    for i = 1:num_methods
        if j <= 2, scaled_metrics(i, j) = (raw_metrics(i, j) - med) / mad_val;
        else, scaled_metrics(i, j) = (med - raw_metrics(i, j)) / mad_val; end
    end
end

ranks_matrix = zeros(num_methods, num_scenarios);
for s = 1:num_scenarios
    w = scenarios{s, 2}; scores = zeros(num_methods, 1);
    for i = 1:num_methods
        scores(i) = sum(w .* scaled_metrics(i,:));
    end
    [~, sort_idx] = sort(scores, 'descend');
    for r = 1:num_methods, ranks_matrix(sort_idx(r), s) = r; end
end

fprintf('=> 数据计算完成！开始绘制绝对排版图表...\n');

% ================= 开始画图 =================
fig = figure('Position', [100, 100, 1050, 520], 'Color', 'w');

% 使用绝对坐标，防止图例挤压画布
ax = axes('Position', [0.1, 0.15, 0.65, 0.75]);
hold(ax, 'on');

% 【定制专属颜色】: 与图1完全一致
base_colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'}; % 蓝(BO), 橙(GA), 黄(PSO), 紫(SSA)
markers_pool = {'o', 's', '^', 'd'};

for m = 1:num_methods
    p_name = valid_methods{m};
    
    % 识别基础算法赋予对应颜色
    if contains(p_name, 'bayes'), color_idx = 1; base_name = 'BO';
    elseif contains(p_name, 'ga'), color_idx = 2; base_name = 'GA';
    elseif contains(p_name, 'pso'), color_idx = 3; base_name = 'PSO';
    elseif contains(p_name, 'ssa'), color_idx = 4; base_name = 'SSA';
    end
    
    % 混合算法用粗实线，独立算法用细虚线
    if contains(p_name, '_f')
        lw = 2.8; ls = '-'; 
        disp_name = [base_name '-Fmincon'];
    else
        lw = 1.6; ls = '--'; 
        disp_name = base_name;
    end
    
    plot(ax, 1:num_scenarios, ranks_matrix(m, :), 'Color', base_colors{color_idx}, 'LineWidth', lw, ...
        'LineStyle', ls, 'Marker', markers_pool{color_idx}, 'MarkerSize', 10, ...
        'MarkerFaceColor', 'w', 'DisplayName', disp_name);
end

% 【核心修正 1：Y 轴强力反转与锁死】第一名绝对在最上方
set(ax, 'YDir', 'reverse');
ylim(ax, [0.5, 8.5]); 
% 强力剥夺 MATLAB 生成副刻度的权限
set(ax, 'YTick', 1:8, 'YTickMode', 'manual');

% 【核心修正 2：X 轴刻度强力锁死】彻底消灭 1.5, 2.5 小数
xlim(ax, [0.5, 5.5]);
set(ax, 'XTick', 1:5, 'XTickMode', 'manual'); 
set(ax, 'XTickLabel', scenario_labels);

% 【字体修改】：全图 15号 新罗马字体
set(ax, 'FontSize', 15, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.2);

ylabel(ax, '$\mathbf{SGPI_{obj}}$ \textbf{rank (1=Best)}', 'Interpreter', 'latex', 'FontSize', 15);
xlabel(ax, 'Weighting Scenarios', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% 【核心修正 3：图注生成】放在画布外部右侧，安全不遮挡
legend(ax, 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 15, 'FontName', 'Times New Roman');

grid(ax, 'on'); ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
hold(ax, 'off');

% 设置保存目录（修改为当前目录下的 Figures 文件夹）
save_dir = fullfile(pwd, 'Figures');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% 拼接完整路径
save_path = fullfile(save_dir, 'Figure2_SGPI_Weight_Sensitivity.emf');

% 使用 -dmeta 参数导出 EMF（MATLAB 2019b 完美兼容，最适合插入 Word）
print(fig, save_path, '-dmeta', '-r600');
fprintf('=> [图3] SGPI权重敏感性图 已成功保存至: %s\n', save_path);