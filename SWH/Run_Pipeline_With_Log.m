clear; clc;
% 1. 强制加载当前目录及其所有子目录（确保所有函数都在路径中）
addpath(genpath(pwd));

% 2. 定义日志文件名（带当前运行时间戳，防止多次运行被覆盖）
log_filename = sprintf('SWH_Pipeline_FullRunLog_%s.txt', datestr(now, 'yyyymmdd_HHMMSS'));

% === 修复核心：将变量存入临时文件避险，防止被 main.m 里的 clear 删掉 ===
save('temp_log_filename_backup.mat', 'log_filename');

% 3. 开启日志记录仪
diary(log_filename);

fprintf('===================================================================\n');
fprintf('  启动全局日志监控：所有命令行输出将同步写入本地文本\n');
fprintf('  日志文件名称：%s\n', log_filename);
fprintf('===================================================================\n\n');

try
    % 3. 运行主优化程序
    if exist('main.m', 'file')
        fprintf('[流程 1/2] 正在启动主优化率定程序 (main.m)...\n');
        fprintf('-------------------------------------------------------------------\n');
        run('main.m'); 
        fprintf('\n-------------------------------------------------------------------\n');
        fprintf('[流程 1/2] 主优化率定程序运行完毕。\n');
    else
        fprintf('[警告] 未在当前目录下找到 main.m 主程序，跳过此步。\n');
    end
    
    % 4. 运行后处理与重修指标脚本（注释已解除）
    if exist('generate_rebuttal_exports_and_metrics.m', 'file')
        fprintf('\n[流程 2/2] 正在启动后处理与统计脚本 (generate_rebuttal_exports_and_metrics.m)...\n');
        fprintf('-------------------------------------------------------------------\n');
        run('generate_rebuttal_exports_and_metrics.m');
        fprintf('\n-------------------------------------------------------------------\n');
        fprintf('[流程 2/2] 后处理与重修指标脚本运行完毕。\n');
    else
        fprintf('[警告] 未找到 generate_rebuttal_exports_and_metrics.m，跳过此步。\n');
    end

catch ME
    % 捕捉运行中可能出现的任何中断或报错
    fprintf('\n===================================================================\n');
    fprintf(' [错误提示] 程序在运行期间发生异常中断：\n');
    fprintf(' 错误信息: %s\n', ME.message);
    fprintf(' 出错位置: 行 %d (%s)\n', ME.stack(1).line, ME.stack(1).name);
    fprintf('===================================================================\n');
end

% 5. 关闭日志记录仪
diary off;

% === 修复核心：所有子脚本运行完毕后，从避难所读取回 log_filename，并销毁临时文件 ===
if exist('temp_log_filename_backup.mat', 'file')
    load('temp_log_filename_backup.mat', 'log_filename');
    delete('temp_log_filename_backup.mat'); % 用完即焚，保持目录干净
end

fprintf('\n===================================================================\n');
fprintf('  全流程自动运行结束！\n');
fprintf('  完整的运行日志（含4算法过程、生育期结果、时间统计、Top10指标）已保存至：\n');
fprintf('  >> %s <<\n', fullfile(pwd, log_filename));
fprintf('===================================================================\n');