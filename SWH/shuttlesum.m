% 修复后的 shuttlesum.m
function [shuttledaily,shuttlemonthly] = shuttlesum(shuttle_halfhour, file_suffix)
 % ========== 新增这一行，防止未赋值报错 ==========
    shuttlemonthly = []; 

    if nargin < 2
        file_suffix = ''; % 默认无后缀
    end
    shuttle = shuttle_halfhour;    
    [m, n] = size(shuttle);
    % --------------------------- daily ---------------------------------------
    daylength = floor(m/48);
    shuttledaily = zeros(daylength, 23); % 修复：扩展到23列以容纳E/ET
    dayl = 30; % the minimum data in a day

    for i = 1:daylength
        shuttledaily(i,1:3) = shuttle((i-1)*48+1, 1:3);
        shuttledaily(i,4)   = shuttle((i-1)*48+1, 5);
        
        dayflux = zeros(48, 22);
        dayflux(:,5:9)   = shuttle((i-1)*48+1:i*48, 7:11); % Ta, RH, VPD, Canopy height, SW
        dayflux(:,10:11) = shuttle((i-1)*48+1:i*48, 14:15); % wind speed, CO2
        dayflux(:,12)    = shuttle((i-1)*48+1:i*48, 19);    % LAI
        dayflux(:,13:14) = shuttle((i-1)*48+1:i*48, 25:26); % rsc, rss
        dayflux(:,15:16) = shuttle((i-1)*48+1:i*48, 12:13); % Rn, G
        dayflux(:,17:19) = shuttle((i-1)*48+1:i*48, 16:18); % GPP, measured ET, Rain 
        dayflux(:,20:22) = shuttle((i-1)*48+1:i*48, 21:23); % T, E, modelled ET 

        % for mean values
        for k = 5:12
            valid_data = dayflux(:,k);
            valid_data = valid_data(valid_data ~= -99999);
            if length(valid_data) < dayl
               sumflux = -99999;     % 修复了 sumfux 的拼写错误
            else
               sumflux = mean(valid_data);         
            end
            shuttledaily(i,k) = sumflux;
        end
        
        % rsc and rss during daytime 10:00-16:00 (rows 21-33)
        for k = 13:14
            valid_data = dayflux(21:33, k);
            valid_data = valid_data(valid_data ~= -99999 & valid_data ~= 0);
            if length(valid_data) < 3
               sumflux = -99999;     
            else
               sumflux = mean(valid_data);         
            end
            shuttledaily(i,k) = sumflux;
        end
           
        % sum values (Rn, G)
        for k = 15:16
            valid_data = dayflux(:,k);
            valid_data = valid_data(valid_data ~= -99999);
            if length(valid_data) < dayl
               sumflux1 = -99999;     
            else
               sumflux1 = mean(valid_data) * 1800 * 48 / 1000000;        
            end
            shuttledaily(i,k) = sumflux1;
        end
        
        % for GPP
        valid_data = dayflux(:,17);
        valid_data = valid_data(valid_data ~= -99999);
        if length(valid_data) < dayl
           sumflux1 = -99999;     
        else
           sumflux1 = mean(valid_data) * 23.56363636;  % mg co2/m2.s to g C/m2.day   
        end
        shuttledaily(i,17) = sumflux1;
            
        % measured ET, Rain, Modeled T, E, ET
        for k = 18:22
            valid_data = dayflux(:,k);
            valid_data = valid_data(valid_data ~= -99999);
            if length(valid_data) < dayl
               sumflux1 = -99999;     
            else
               if k == 19 % Rain
                   sumflux1 = mean(valid_data) * 48; 
               else % ET components
                   sumflux1 = mean(valid_data) * 1.8 * 48;
               end
            end
            shuttledaily(i,k) = sumflux1;
        end
    end
    
    % E/ET Ratio
    ratio = shuttledaily(:,21) ./ shuttledaily(:,22);
    ratio(isnan(ratio) | isinf(ratio)) = -99999; % 防止除以0报错
    shuttledaily(:,23) = ratio;

    % 写入Excel (添加 .xlsx 后缀防止COM错误)
    headerlineDaily = {'Year' 'Month' 'Day' 'DOY' 'Ta' 'RH' 'VPD' 'Canopy height' 'SW' 'wind speed' 'CO2 concentration' 'LAI' 'rsc' 'rss' 'Rn' 'G' 'GPP' 'measured ET' 'Rain' 'Modeled T' 'Modeled E' 'Modeled ET' 'E/ET'}; 
    filename = 'shuttledaily';
    if ~isempty(file_suffix), filename = [filename '_' file_suffix]; end
    try
        xlswrite([filename '.xlsx'], headerlineDaily, 'sheet1', 'A1'); 
        xlswrite([filename '.xlsx'], shuttledaily, 'sheet1', 'A2');
    catch
        disp('Excel写入失败，已跳过。');
    end


end