%----------------------------------------------------------
%Rn,G,vpd,RH,Ta, Ca,U,LAI,GPP,SW,windspeed
% estimate of rss, rsc have great uncertainty,
% make sure the units consistent!!!!!!!!!, special attention to rsc
% Calculation of rss AND rsc is SPECIAL at SD!!!!!
% ---------------------------------------------------------

 function [shuttle_halfhour]=halfhour_shuttleworth_fillgap(shuttle_FlagBadData,para_range_10,Qs)
 
 shuttle=shuttle_FlagBadData;
 %ma=1; 注意不要与m=length(shuttle)的m混淆

  B1=para_range_10(3,1); %均值位于第三行
  B2=para_range_10(3,2);
  B3=para_range_10(3,3); % B1, B2,B3 are the coeficient for estimating rss, Mahfouf's method. unite of B3: m s-1
  a1=para_range_10(3,4); % go,a1 is  the coeficient for estimating rsc, Ball-Berry's method,
  g01=para_range_10(3,5);  % mol m-2 s-1, not mol m-2 s-1
  LightExtCoef=para_range_10(3,6);  % 
  
    h=shuttle(:,10); % cnopy height 0.4 m
   
   %  reference height x
x = zeros(size(h));  % 创建与 h 同大小的列向量
w = zeros(size(h));  % 创建与 h 同大小的列向量

o = find(h <= 1);
x(o) = 2;
w(o) = h(o) ./ 10;
clear o;
o = find(h > 1 & h <= 3);
x(o) = 1.5 + h(o);
w(o) = h(o) ./ 20;
clear o;
o = find(h > 3);
x(o) = 1.5 * h(o);
w(o) = 0.09;%-------------
clear o ;

  z01=0.01;%  effective lenght=0.01 m
  d=h*0.63; % zero displacement
  z0=h*0.13; % roughness
  k=0.41;

  Cp=1017.388; % [J/(kg K)].
  d_a=1.1723;  % g/L=kg.m-3
  r=0.6721;      % psychrometric constant mb K-1 or hPa K-1/
  
  %--------------------
  Rn=shuttle(:,12); % w s-1
  G=shuttle(:,13);  % w s-1
  VPD=shuttle(:,9)*10; % 10 is the coversion coefficient from kPa to hPa
  RH=shuttle(:,8)/100; 
  Ta=shuttle(:,7);
  Ca=shuttle(:,15)*0.509; % mg Co2 m-3 to umol mol-1
  U=shuttle(:,14)+0.0001; % wind speed, must not be zero when run the PEST
  LAI=shuttle(:,19)+0.0001; % must not be zero when run the PEST
  GPP=shuttle(:,16);
  SW=shuttle(:,11); % soil moisture at 5 cm
[slp]=(4098./(Ta+237.3).^2).*(0.6108*exp(17.27*Ta./(237.3+Ta)))*10; % 10 is the conversion coeficient of slp from kPa/K to hPa/K
% save slp.dat slp -ascii;
 % --------soil and canopy surface resistance--------- 
 %rss=B3;
 [rss]=B1*(Qs./SW).^B2+B3; %+2.5*litter; % litter is the dead matter on the bare soil;Mahfouf J and Noihan J.1991.comparative study of various formulations of evaporation from bare soil using in situ data. Journal of applied meteorology, 30:1354-1365.
 % values of in the original paper are B1=3.5, B2=2.3,B3=33.5
% rsc=800;
 g0=g01*LAI;
[gsw]=(g0+a1*22.73*GPP.*RH./Ca)*0.0224*1.6; %22.73 is the conversion coefficient of GPP from mg CO2 m-2 s-1 to umol m-2 s-1, 0.0224 is the conversion coefficient of canopy conductance from mol m-2 s-1 to m s-1                                                                                  
%[gsw]=(g0+a1*22.73*GPP.*RH./Ca)*0.0224*1.6;
[rsc]=(1./gsw); % wang,Ball-Berry 
 
% fQ=1-exp(-Q/500);
% fVPD=1-fcvpd*VPD;
% fT=1-fct*(25-Ta).^2;
% rsc=((0.5*LAI+1)./LAI)*r0./(fQ.*fVPD.*fT); % Olioso et al., 1996

%  save rsc.dat rsc -ascii; 

% ---------arerodynamic resistance--------- 
[ras0]=(log(x./z01)).*(log((d+z0)./z01)./(k^2*U));
[raa0]=((log(x./z01)).^2)./(k^2*U)-ras0;
  % z01: effective roughness length of the substrate (0.01 m).
[raaa]=log((x-d)./z0)./(k^2*U).*(log((x-d)./(h-d))+h./(2.5*(h-d)).*(exp(2.5*(1-(d+z0)./h))-1)); 
[rasa]=(log((x-d)./z0))./(k^2*U).*(h./(2.5*(h-d))).*(exp(2.5)-exp(2.5*(1-(d+z0)./h)));
  % x: reference height(2 m),d:zero displacement(0.63h),z0:canopy roughness
  % length (0.13h); k: karman constant (0.41); U: wind speed at the height
  % of 2 m; h: height of canopy; 
% 对于稀疏冠层（0<LAI<4）有：
% 修复条件判断
index = LAI < 4;
raa = zeros(size(LAI));
ras = zeros(size(LAI));

raa(index) = 0.25*LAI(index).*raaa(index) + 0.25*(4-LAI(index)).*raa0(index);
ras(index) = 0.25*LAI(index).*rasa(index) + 0.25*(4-LAI(index)).*ras0(index);

raa(~index) = raaa(~index);
ras(~index) = rasa(~index);
% save raa.dat raa -ascii; 
% save ras.dat ras -ascii; 

% ---------rac----------------------
Rns=Rn.*exp(-LightExtCoef*LAI);
rb=130*(w./U).^0.5; % w: width of leaf; U: wind speed at the top of canopy  in Yu's review paper, choose Monteith's number
rac=rb./(2*LAI);


[Ra]=(slp+r).*raa;
[Rs]=(slp+r).*ras+r.*rss;
[Rc]=(slp+r).*rac+r.*rsc;

Cc=1./(1+Rc.*Ra./(Rs.*(Rc+Ra)));
Cs=1./(1+Rs.*Ra./(Rc.*(Rs+Ra)));

% save Ra.dat Ra -ascii;
% save Rs.dat Rs -ascii;
% save Rc.dat Rc -ascii;
% save Cc.dat Cc -ascii;
% save Cs.dat Cs -ascii;

[fenzi1]=slp.*(Rn-G)+(d_a*Cp*VPD-slp.*rac.*(Rns-G))./(raa+rac);
[fenmu1]=slp+r*(1+rsc./(raa+rac));

[PMc]=fenzi1./fenmu1;

[fenzi2]=slp.*(Rn-G)+(d_a*Cp*VPD-slp.*ras.*(Rn-Rns))./(raa+ras);
[fenmu2]=slp+r*(1+rss./(raa+ras));
[PMs]=fenzi2./fenmu2;

[Ec]=Cc.*PMc/2450; % 1wa/s=1 J, 1 g H2O=2450 J
[Es]=Cs.*PMs/2450;
[E]=Ec+Es;
[Es_E]=Es./E;

shuttle(:,21)=Ec;
shuttle(:,22)=Es;
shuttle(:,23)=E;
shuttle(:,24)=Es_E;
shuttle(:,25)=rsc;
shuttle(:,26)=rss;

  clear  Rn;
  clear G;
  clear VPD;
  clear RH;
  clear Ta;
  clear Ca; % co2 ppmv
  clear U; % wind speed
  clear LAI;
  clear GPP;
  clear SW;
  clear slp;
  

  %--------------------QA/QC ----------------------------------------
  m=length(shuttle(:,1));
    for j=1:m
        for k=7:17
          if shuttle(j,20)==k
            shuttle(j,k)=-99999;
            shuttle(j,21:26)=-99999;
          end
        end 
    end    
 %--------------------Gap filling----------------------------------------
  
      shuttle_halfhour=shuttle(:,1:26); % shuttle_halfhour will be the data with input variables being gap filled
     MDVinsert=zeros(m,11);
     MDVinsert(:,1:11)=shuttle_halfhour(:,7:17); % Ta 	RH 	VPD 	PPFD	SW_5cm	Rn	G 	windspeed	co2 	GPP	 ET_measure E_measure


%________linear filling___________

for k=1:11
    for i=2:m
        if MDVinsert(i,k)==-99999&MDVinsert(i-1,k)~=-99999&MDVinsert(i+1,k)~=-99999
            MDVinsert(i,k)=(MDVinsert(i-1,k)+MDVinsert(i+1,k))/2;
        end
    end
end

%___________MDV filling______________
%For the medial days  
for k=1:11
    for i=3*48+1:m-3*48
        if MDVinsert(i,k)==-99999
           for j=-3:3
               tempcrop(j+4,1)=MDVinsert(i+j*48,k);
           end
           o=tempcrop~=-99999;
           io=find(o);
           if length(tempcrop(io))>2
              MDVinsert(i,k)=mean(tempcrop(io));
           end
        end
    end
%For the first 3 days 
    for i=1:3*48
        if MDVinsert(i,k)==-99999
           for j=1:7
               tempcrop(j,1)=MDVinsert(i+j*48,k);
           end
           o=tempcrop~=-99999;
           io=find(o);
           if length(tempcrop(io))>2
              MDVinsert(i,k)=mean(tempcrop(io));
           end
        end
    end
% For the last 3 days
    for i=m-3*48:m
        if MDVinsert(i,k)==-99999
           for j=1:7
               tempcrop(j,1)=MDVinsert(i-j*48,k);
           end
           o=tempcrop~=-99999;
           io=find(o);
           if length(tempcrop(io))>2
              MDVinsert(i,k)=mean(tempcrop(io));
           end
        end
    end
end
%

shuttle_halfhour(:,7:17)=MDVinsert(:,1:11);            

%

clear k;
clear i;
clear j;
clear MDVinsert;




%______________re_estimate E, T and ET____________________________________
%   
  B1=para_range_10(3,1); %均值位于第三行
  B2=para_range_10(3,2);
  B3=para_range_10(3,3); % B1, B2,B3 are the coeficient for estimating rss, Mahfouf's method. unite of B3: m s-1
  a1=para_range_10(3,4);  % go,a1 is  the coeficient for estimating rsc, Ball-Berry's method,
  g01=para_range_10(3,5);  % mol m-2 s-1, not mol m-2 s-1
  LightExtCoef=para_range_10(3,6);
  
   h=shuttle_halfhour(:,10); % cnopy height 0.4 m
   
   %  reference height x
x = zeros(size(h));  % 创建与 h 同大小的列向量
w = zeros(size(h));  % 创建与 h 同大小的列向量

o = find(h <= 1);
x(o) = 2;
w(o) = h(o) ./ 10;
clear o;
o = find(h > 1 & h <= 3);
x(o) = 1.5 + h(o);
w(o) = h(o) ./ 20;
clear o;
o = find(h > 3);
x(o) = 1.5 * h(o);
w(o) = 0.09;%-------------
clear o ;
  z01=0.01;%  effective lenght=0.01 m
  d=h*0.63; % zero displacement
  z0=h*0.13; % roughness
  k=0.41;

  Cp=1017.388; % [J/(kg K)].
  d_a=1.1723;  % g/L=kg.m-3
  r=0.6721;      % psychrometric constant mb K-1 or hPa K-1/

  %--------------------
  Rn=shuttle_halfhour(:,12); % w s-1
  G=shuttle_halfhour(:,13);  % w s-1
  VPD=shuttle_halfhour(:,9)*10; % 10 is the coversion coefficient from kPa to hPa
  RH=shuttle_halfhour(:,8)/100; 
  Ta=shuttle_halfhour(:,7);
  Ca=shuttle_halfhour(:,15)*0.509; % mg Co2 m-3 to umol mol-1
  U=shuttle_halfhour(:,14)+0.0001; % wind speed, must not be zero when run the PEST
  LAI=shuttle_halfhour(:,19)+0.0001; % must not be zero when run the PEST
  GPP=shuttle_halfhour(:,16);
  SW=shuttle_halfhour(:,11); % soil moisture at 5 cm
[slp]=(4098./(Ta+237.3).^2).*(0.6108*exp(17.27*Ta./(237.3+Ta)))*10; % 10 is the conversion coeficient of slp from kPa/K to hPa/K
% save slp.dat slp -ascii;
 % --------soil and canopy surface resistance--------- 
% rss=1000;
%[rss]=A*exp(B*SW)+2.5*litter
%[rss]=B3;
[rss]=B1*(Qs./SW).^B2+B3; %+2.5*litter; % litter is the dead matter on the bare soil;Mahfouf J and Noihan J.1991.comparative study of various formulations of evaporation from bare soil using in situ data. Journal of applied meteorology, 30:1354-1365.
 % values of in the original paper are B1=3.5, B2=2.3,B3=33.5
% rsc=800;
g0=g01*LAI;
[gsw]=(g0+a1*22.73*GPP.*RH./Ca)*0.0224*1.6;
%[gsw]=(g0+a1*22.73*GPP.*RH.*fSW./Ca)*0.0224; %22.73 is the conversion coefficient of GPP from mg CO2 m-2 s-1 to umol m-2 s-1, 0.0224 is the conversion coefficient of canopy conductance from mol m-2 s-1 to m s-1                                                                                  
[rsc]=(1./gsw); % wang,Ball-Berry 
 
% fQ=1-exp(-Q/500);
% fVPD=1-fcvpd*VPD;
% fT=1-fct*(25-Ta).^2;
% rsc=((0.5*LAI+1)./LAI)*r0./(fQ.*fVPD.*fT); % Olioso et al., 1996

%  save rsc.dat rsc -ascii; 

% ---------arerodynamic resistance--------- 
[ras0]=(log(x./z01)).*(log((d+z0)./z01)./(k^2*U));
[raa0]=((log(x./z01)).^2)./(k^2*U)-ras0;
  % z01: effective roughness length of the substrate (0.01 m).
[raaa]=log((x-d)./z0)./(k^2*U).*(log((x-d)./(h-d))+h./(2.5*(h-d)).*(exp(2.5*(1-(d+z0)./h))-1)); 
[rasa]=(log((x-d)./z0))./(k^2*U).*(h./(2.5*(h-d))).*(exp(2.5)-exp(2.5*(1-(d+z0)./h)));
  % x: reference height(2 m),d:zero displacement(0.63h),z0:canopy roughness
  % length (0.13h); k: karman constant (0.41); U: wind speed at the height
  % of 2 m; h: height of canopy; 
% 对于稀疏冠层（0<LAI<4）有：
% 修复条件判断
index = LAI < 4;
raa = zeros(size(LAI));
ras = zeros(size(LAI));

raa(index) = 0.25*LAI(index).*raaa(index) + 0.25*(4-LAI(index)).*raa0(index);
ras(index) = 0.25*LAI(index).*rasa(index) + 0.25*(4-LAI(index)).*ras0(index);

raa(~index) = raaa(~index);
ras(~index) = rasa(~index);
% save raa.dat raa -ascii; 
% save ras.dat ras -ascii; 

% ---------rac----------------------
Rns=Rn.*exp(-LightExtCoef*LAI);
rb=130*(w./U).^0.5; % w: width of leaf; U: wind speed at the top of canopy  in Yu's review paper, choose Monteith's number
rac=rb./(2*LAI);


[Ra]=(slp+r).*raa;
[Rs]=(slp+r).*ras+r.*rss;
[Rc]=(slp+r).*rac+r.*rsc;

Cc=1./(1+Rc.*Ra./(Rs.*(Rc+Ra)));
Cs=1./(1+Rs.*Ra./(Rc.*(Rs+Ra)));

% save Ra.dat Ra -ascii;
% save Rs.dat Rs -ascii;
% save Rc.dat Rc -ascii;
% save Cc.dat Cc -ascii;
% save Cs.dat Cs -ascii;

[fenzi1]=slp.*(Rn-G)+(d_a*Cp*VPD-slp.*rac.*(Rns-G))./(raa+rac);
[fenmu1]=slp+r*(1+rsc./(raa+rac));

[PMc]=fenzi1./fenmu1;

[fenzi2]=slp.*(Rn-G)+(d_a*Cp*VPD-slp.*ras.*(Rn-Rns))./(raa+ras);
[fenmu2]=slp+r*(1+rss./(raa+ras));
[PMs]=fenzi2./fenmu2;

[Ec]=Cc.*PMc/2450; % 1wa/s=1 J, 1 g H2O=2450 J
[Es]=Cs.*PMs/2450;
[E]=Ec+Es;
[Es_E]=Es./E;

shuttle_halfhour(:,21)=Ec;
shuttle_halfhour(:,22)=Es;
shuttle_halfhour(:,23)=E;
shuttle_halfhour(:,24)=Es_E;
shuttle_halfhour(:,25)=rsc;
shuttle_halfhour(:,26)=rss;

  clear  Rn;
  clear G;
  clear VPD;
  clear RH;
  clear Ta;
  clear Ca; % co2 ppmv
  clear U; % wind speed
  clear LAI;
  clear GPP;
  clear SW;
  clear slp;
  
  %---------gap fill again-----------------
 MDVinsert(:,1:3)=shuttle_halfhour(:,21:23); % T,E,Emodel
     

          %________linear filling___________

for k=1:3
    for i=2:m
        if MDVinsert(i,k)==-99999&MDVinsert(i-1,k)~=-99999&MDVinsert(i+1,k)~=-99999
            MDVinsert(i,k)=(MDVinsert(i-1,k)+MDVinsert(i+1,k))/2;
        end
    end
end

         %___________MDV filling______________
for k=1:3
  %for media days 
    for i=3*48+1:m-3*48
        if MDVinsert(i,k)==-99999
           for j=-3:3
               tempcrop(j+4,1)=MDVinsert(i+j*48,k);
           end
           o=tempcrop~=-99999;
           io=find(o);
           if length(tempcrop(io))>2
              MDVinsert(i,k)=mean(tempcrop(io));
           end
        end
    end
%For the first 3 days    
    for i=1:3*48
        if MDVinsert(i,k)==-99999
           for j=1:7
               tempcrop(j,1)=MDVinsert(i+j*48,k);
           end
           o=tempcrop~=-99999;
           io=find(o);
           if length(tempcrop(io))>2
              MDVinsert(i,k)=mean(tempcrop(io));
           end
        end
    end
% For the last 3 days
    for i=m-3*48:m
        if MDVinsert(i,k)==-99999
           for j=1:7
               tempcrop(j,1)=MDVinsert(i-j*48,k);
           end
           o=tempcrop~=-99999;
           io=find(o);
           if length(tempcrop(io))>2
              MDVinsert(i,k)=mean(tempcrop(io));
           end
        end
    end
end
%

shuttle_halfhour(:,21:23)=MDVinsert(:,1:3);
shuttle_halfhour(:,24)=shuttle_halfhour(:,22)./shuttle_halfhour(:,23);

   
headerlineHourly={'Year' 'Month' 'Day' 'Hour' 'DOY' 'Don' 'Ta' 'RH' 'VPD' 'Canopy height' 'SW' 'Rn' 'G' 'wind speed' 'CO2 concentration' 'GPP'...
     'measured ET' 'Rain' 'LAI'  'DataFlag'  'Modeled T' 'Modeled E' 'Modeled ET' 'E/ET' 'rsc' 'rss'}; 
 xlswrite('shuttle_halfhour',headerlineHourly,'sheet1','A1'); 
 xlswrite('shuttle_halfhour',shuttle_halfhour,'sheet1','A2');
    % save shuttle_halfhour.dat shuttle_halfhour -ascii;
     clear shuttle;

