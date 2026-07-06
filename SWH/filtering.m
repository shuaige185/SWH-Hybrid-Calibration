%---------------------------------------------------------------------
% The program is used to filter the bad-quality data before estimating parameters
%       1 Precipitation
%       2 CO2 concentration 
%       3 Threshold
%       4 3 times standard deviation

%---------------------------------------------------------------------
  %SW2=dlmread('E:\SW.dat', '\t');  
  function [shuttle_backup,shuttle_FlagBadData]=filtering(shuttle_original)
  SW2=shuttle_original;
   
    m=length(SW2(:,1));
    for i=1:m    
        if SW2(i,7)>60|SW2(i,7)<-50         % Ta
           SW2(i,20)=7;
        end  
        
       if SW2(i,8)>100|SW2(i,8)<0         % RH
           SW2(i,20)=8;
       end   
        
        if SW2(i,9)<0|SW2(i,9)>10   %vpd
            SW2(i,20)=9;
        end
        
        if SW2(i,10)>100|SW2(i,10)<0         % Canopy height
           SW2(i,20)=10;
        end   
        
        if SW2(i,11)>0.7|SW2(i,11)<=0         % Soil water content
           SW2(i,20)=11;
        end   
        
        if SW2(i,12)<-200|SW2(i,12)>1000 %Rn   
            SW2(i,20)=12;
        end
         
        if SW2(i,13)<-100|SW2(i,13)>400 %G   
            SW2(i,20)=13;
        end
        
        if SW2(i,14)<0|SW2(i,14)>100 %Wind speed
            SW2(i,20)=14;
        end  
        
        if SW2(i,15)>1200|SW2(i,15)<300         % CO2 concentration     
           SW2(i,20)=15;
        end   
        
        if SW2(i,16)<-1|SW2(i,16)>3 %GPP site-specific threshold, duolun -0.2~0.3
            SW2(i,20)=16;
        end

        if SW2(i,17)<-1|SW2(i,17)>4 %ET gap-filled by MDV or ET-Rn relationship, duolun -0.1~0.45
            SW2(i,20)=17;
        end

%         if SW2(i,19)>0&SW2(i,19)<100;   % if there is precipitation the flux is labeled as -99999
%            SW2(i,21)=19;                        %
%         end;  
        
%         if SW2(i,5)<=160&SW2(i,7)<-0.4&SW2(i,7)~=-99999; % mark the abnormal daytime flux in May and June
%            SW2(i,7)=-99999;
%        end;
%         if SW2(i,5)>=275&SW2(i,7)<-0.4&SW2(i,7)~=-99999; % mark the abnormal daytime flux in September
%            SW2(i,7)=-99999;
%        end;
      
    end
    
 shuttle_FlagBadData=SW2; % 
% 
% save shuttle_flagbaddata.dat shuttle_flagbaddata -ascii;

%-------------筛选合格数据用于模型参数估计---------------
    o=SW2(:,20)==0;
    io=o;
    shuttle_backup=SW2(io,:); 
%     shuttle_carlibration(:,1)=SW2(1:length(shuttle_carlibration),1); % 筛选时可能将输入参数剔掉。

%save shuttle_carlibration.dat shuttle_carlibration -ascii;
