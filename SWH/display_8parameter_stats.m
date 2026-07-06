function display_8parameter_stats(all_method_params, shuttle_calibration, Qs)
% 显示8组参数的统计指标（完美融合MAE并优化排版）
    
    fprintf('\n========== 8组参数统计指标对比 (校准集) ==========\n');
    % 重新优化字符对齐间距，加入 MAE 字段
    fprintf('%-15s %-8s %-8s %-8s %-8s %-8s\n', '参数集', '斜率', 'R?', 'NSE', 'RMSE', 'MAE');
    fprintf('%-15s %-8s %-8s %-8s %-8s %-8s\n', '-----', '----', '---', '---', '----', '---');
    
    % 方法列表
    method_list = {'bayes', 'ga', 'pso', 'ssa'};
    
    for m_idx = 1:length(method_list)
        method = method_list{m_idx};
        
        % 检查两种参数集（原始参数与混合优化参数）
        for suffix_idx = 1:2
            if suffix_idx == 1
                suffix = '';
                param_name = method;
            else
                suffix = '_f';
                param_name = [method suffix];
            end
            
            if isfield(all_method_params, param_name) && ...
               ~isempty(all_method_params.(param_name))
                
                params = all_method_params.(param_name);
                try
                    % 直接调用全局核心模型进行半小时尺度模拟
                    [shuttle_temp] = SWH_halfhour(shuttle_calibration, params(1), params(2), params(3), params(4), params(5), params(6), Qs);
                    ET_sim = shuttle_temp(:,23);
                    ET_obs = shuttle_calibration(:,17);
                    
                    % 严格过滤无效观测数据
                    valid_idx = ET_obs ~= -99999 & ~isnan(ET_sim) & ~isnan(ET_obs);
                    ET_sim = ET_sim(valid_idx);
                    ET_obs = ET_obs(valid_idx);
                    
                    if length(ET_sim) >= 50
                        % 1. 计算斜率 (Slope)
                        p_ET = polyfit(ET_obs, ET_sim, 1);
                        slope = p_ET(1);
                        
                        % 2. 计算决定系数 (R?)
                        coeCorr = corrcoef(ET_sim, ET_obs);
                        r2 = coeCorr(2,1)^2;
                        
                        % 3. 计算均方根误差 (RMSE)
                        rmse = sqrt(mean((ET_sim - ET_obs).^2));
                        
                        % 4. 计算平均绝对误差 (MAE)
                        mae = mean(abs(ET_sim - ET_obs));
                        
                        % 5. 计算纳什效率系数 (NSE)
                        numerator = sum((ET_sim - ET_obs).^2);
                        denominator = sum((ET_obs - mean(ET_obs)).^2);
                        if denominator == 0, nse = -Inf; else, nse = 1 - (numerator / denominator); end
                    else
                        slope = NaN; r2 = NaN; nse = NaN; rmse = NaN; mae = NaN;
                    end
                catch
                    slope = NaN; r2 = NaN; nse = NaN; rmse = NaN; mae = NaN;
                end
                
                % 严格按照上方表头的间距控制，格式化输出浮点数
                fprintf('%-15s %-8.3f %-8.3f %-8.3f %-8.4f %-8.4f\n', ...
                    param_name, slope, r2, nse, rmse, mae);
            end
        end
    end
    fprintf('==================================================\n');
end