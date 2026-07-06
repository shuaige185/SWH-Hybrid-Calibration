% 新建文件：algorithm_config.m
function config = algorithm_config(algorithm_name, data_size)
    % 为不同算法提供优化配置
    
    config = struct();
    
    switch lower(algorithm_name)
        case 'ga'
            % 遗传算法配置
            config.population_size = 200;
            config.max_iterations = 250;
            config.elite_count = 5;
            config.crossover_rate_initial = 0.9;
            config.mutation_rate_initial = 0.3;
            config.selection_pressure = 1.5;
            config.sbx_eta = 20;
            config.pm_eta = 20;
            
        case 'pso'
            % 粒子群算法配置
            config.population_size = 200;
            config.max_iterations = 250;
            config.w_max = 0.9;
            config.w_min = 0.4;
            config.c1_max = 2.5;
            config.c1_min = 1.5;
            config.c2_max = 2.5;
            config.c2_min = 1.5;
            config.topology = 'ring';
            config.neighborhood_size = 5;
            config.v_max_factor = 0.15;
            
        case 'ssa'
            % 麻雀算法配置
            config.population_size = 200;
            config.max_iterations = 250;
            config.pd = 0.7;        % 发现者比例
            config.sd = 0.2;        % 警戒者比例
            config.ST = 0.6;        % 安全阈值
            config.R2_min = 0.8;    % 预警值最小值
            config.R2_max = 1.0;    % 预警值最大值
            config.alpha = 1.0;     % 发现者更新参数
            
        case 'bayes'
            % 贝叶斯优化配置
            config.max_evaluations = 200;
            config.acquisition_function = 'expected-improvement-plus';
            config.exploration_ratio = 0.6;
            config.num_seed_points = 15;
            
        otherwise
            error('未知算法: %s', algorithm_name);
    end
    
    % 根据数据量调整参数
    if data_size > 3000
        % 针对3504组数据
        config.population_size = max(config.population_size, 200);
        config.max_iterations = max(config.max_iterations, 250);
    elseif data_size > 1000
        config.population_size = max(config.population_size, 100);
        config.max_iterations = max(config.max_iterations, 100);
    else
        config.population_size = max(config.population_size, 50);
        config.max_iterations = max(config.max_iterations, 50);
    end
    
    fprintf('算法配置: %s (数据量: %d)\n', algorithm_name, data_size);
    fprintf('  种群大小: %d\n', config.population_size);
    fprintf('  最大迭代次数: %d\n', config.max_iterations);
end