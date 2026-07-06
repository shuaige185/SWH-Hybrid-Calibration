% =========================================================================
% 策略B 专属提取、独立重算 GPI 与重排脚本 (Strict Strategy B Isolation)
% 作用：提取 Strategy B 数据，并在 8 种算法的孤立组内重新计算 GPI 和排名
% =========================================================================
clear; clc;

file_name = 'Summary_Metrics_Full_Report.xlsx';
if ~exist(file_name, 'file')
    file_name = 'Summary_Metrics_Full_Report.csv'; 
end
if ~exist(file_name, 'file')
    error('未找到总表，请确认文件在当前目录下！');
end

fprintf('正在读取总表数据...\n');
T = readtable(file_name);

% 1. 提取 Strategy = 'B' 的数据行 (8个算法)
is_B = strcmp(strtrim(string(T.Strategy)), 'B');
T_B = T(is_B, :);

if isempty(T_B)
    error('未找到 Strategy B 的数据！');
end

fprintf('成功提取 %d 行策略 B 数据，正在进行 8 算法内部的孤立重算...\n', height(T_B));

% 2. 遍历每个年份和时间尺度，重新独立计算 GPI 和 1-8 排名
years = unique(T_B.Year);
scales = unique(string(T_B.Scale));

for y = 1:length(years)
    for s = 1:length(scales)
        % 定位到当前的 8 个算法
        idx = find((T_B.Year == years(y)) & strcmp(string(T_B.Scale), scales(s)));
        if isempty(idx), continue; end
        
        subT = T_B(idx, :);
        
        % ==========================================================
        % 核心动作：在这 8 个算法内部重新进行无量纲化和中位数计算
        % ==========================================================
        scale_minmax = @(x) (x - min(x)) ./ (max(x) - min(x) + eps); % 归一化到 0-1
        
        % 仅使用这 8 个数据的极值进行缩放 (即 y_ij)
        y_R2   = scale_minmax(subT.R2);
        y_NSE  = scale_minmax(subT.NSE);
        y_RMSE = scale_minmax(subT.RMSE);
        y_MAE  = scale_minmax(subT.MAE);
        
        % 仅计算这 8 个数据的中位数 (即 tilde{y}_j)
        g_R2   = median(y_R2);
        g_NSE  = median(y_NSE);
        g_RMSE = median(y_RMSE);
        g_MAE  = median(y_MAE);
        
        % 严格按照公式计算新 GPI
        % 公式: sum( alpha_j * (median_j - y_ij) )
        % R2/NSE 越大越好(alpha=-1)；RMSE/MAE 越小越好(alpha=1)
        new_GPI = -1 * (g_R2 - y_R2) ...
                + -1 * (g_NSE - y_NSE) ...
                +  1 * (g_RMSE - y_RMSE) ...
                +  1 * (g_MAE - y_MAE);
                
        % 将纯净版的新 GPI 覆写回表格
        T_B.GPI_Score(idx) = new_GPI;
        
        % ==========================================================
        % 重新分配 1 到 8 名的组内名次
        % ==========================================================
        % 正向指标 (降序)
        [~, p] = sort(subT.R2, 'descend'); r=zeros(size(p)); r(p)=1:length(p); T_B.Rank_R2(idx) = r;
        [~, p] = sort(subT.NSE, 'descend'); r=zeros(size(p)); r(p)=1:length(p); T_B.Rank_NSE(idx) = r;
        [~, p] = sort(subT.KGE, 'descend'); r=zeros(size(p)); r(p)=1:length(p); T_B.Rank_KGE(idx) = r;
        [~, p] = sort(new_GPI, 'descend'); r=zeros(size(p)); r(p)=1:length(p); T_B.Rank_GPI(idx) = r;
        
        % 负向指标 (升序)
        [~, p] = sort(subT.RMSE, 'ascend'); r=zeros(size(p)); r(p)=1:length(p); T_B.Rank_RMSE(idx) = r;
        [~, p] = sort(subT.MAE, 'ascend'); r=zeros(size(p)); r(p)=1:length(p); T_B.Rank_MAE(idx) = r;
        [~, p] = sort(abs(subT.Bias), 'ascend'); r=zeros(size(p)); r(p)=1:length(p); T_B.Rank_absBias(idx) = r;
    end
end

% 3. 导出终极纯净版数据表
output_file = 'Table1_Strategy_B_Recalculated_GPI.xlsx';
try
    writetable(T_B, output_file);
    fprintf('\n=======================================================\n');
    fprintf(' => 成功！策略B 的终极报表已生成: %s\n', output_file);
    fprintf(' => GPI 已完全基于 8 种算法的独立中位数 (tilde{y}) 重新计算！\n');
    fprintf(' => 排名已更新为 1 到 8 名。您可以安全地将其用于表 1。\n');
    fprintf('=======================================================\n');
catch
    output_csv = 'Table1_Strategy_B_Recalculated_GPI.csv';
    writetable(T_B, output_csv);
    fprintf('\n=> 成功！已生成: %s\n', output_csv);
end