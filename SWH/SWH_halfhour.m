% S-W模型，读入的驱动变量为一列数据，而不是一个数

% 2012.4.15 gs乘以1.6
% In Ball-Berry model, the unit of gs is mol m-2 s-1 and the unit of An is umol m-2 s-1 The transfermation coefficient from mgCO2 to umol is 22.73
% For light use efficiency, the transfermation coefficient from  g C/mol IPAR to mg CO2/umol IPAR is 0.0037, and the reverse is 270.27, from mg CO2/umol to mol/mol is 22.73 

% 1 year    2 month    3 day    4 hour    5 doy (day of year)    6 don 
% 7 Ta (air temperature, ℃)    8 RH(relative humidity, %)      9 VPD(kpa)
% 10 canopy height(m)           11 SW (soil volumum water content m3m-3)
% 12 Rn (net raidation, Wa m-2)     13 G (soil heat flux, Wa m-2)    14 windspeed (m s-1) 
% 15 co2( mg Co2 m-3)          16 GPP(mg CO2 m-2s-1)         17 ET_measure (gH2O m-2s-1)
% 18 Rain (mm)                 19 LAI             20 DataFlag (set it as 0)

function [shuttle]=SWH_halfhour(shuttle,b1,b2,b3,a1,g01,LightExtCoef,Qs)

h=shuttle(:,10); % cnopy height 0.5 m
Ta=shuttle(:,7);
RH=shuttle(:,8)/100;
VPD=shuttle(:,9)*10; %from kpa to hPa
SW=shuttle(:,11); % soil moisture at 5 cm
Rn=shuttle(:,12);%  Wa s-1
G=shuttle(:,13); %  Wa s-1
Ca=shuttle(:,15)*0.509; % mg Co2 m-3 to umol mol-1 
U=shuttle(:,14)+0.0001; % wind speed, must not be zero
LAI=shuttle(:,19)+0.0001;% must not be zero
GPP=shuttle(:,16);
[slp]=(4098./(Ta+237.3).^2).*(0.6108*exp(17.27*Ta./(237.3+Ta)))*10; % hpa/k


% ---------estimate of reference the width of leaf w
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
w(o) = 0.09;
clear o;
%-------------

z01=0.01;%  effective lenght=0.01 m
d=h*0.63; % zero displacement
z0=h*0.13; % roughness
k=0.41;

Cp=1017.388;   % [J/(kg K)].
d_a=1.1723;  % g/L=kg.m-3
r=0.6721;      % psychrometric constant mb K-1 or hPa K-1/

    %----shuttleworth model

    % save slp.dat slp -ascii;
    % --------soil and canopy surface resistance---------
    % rss=1000;
    %[rss]=A*exp(B*SW)+2.5*litter

    [rss]=b1.*(Qs./SW).^b2+b3; %+2.5*litter; % litter is the dead matter on the bare soil;Mahfouf J and Noihan J.1991.comparative study of various formulations of evaporation from bare soil using in situ data. Journal of applied meteorology, 30:1354-1365.
    % values of in the original paper are B1=3.5,B2=2.3,B3=33.5
    % rsc=800;
    g0=g01*LAI;
    [gsw]=(g0+a1*22.73*GPP.*RH./Ca)*0.0224*1.6; % before using this model, be sure the unite of each term is consistent with that in references! mg/umol=22.73, 
    %[gsw]=(g0+a1*22.73*GPP.*RH.*fSW./Ca)*0.0224*1.6; % 1.6 from mol Co2 m-2s-2 to mol H2O m-2 s-1, 22.73：from mg CO2 to umol CO2,
    [rsc]=(1./gsw); % wang,Ball-Berry

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

    %
    [Ra]=(slp+r).*raa;
    [Rs]=(slp+r).*ras+r*rss;
    [Rc]=(slp+r).*rac+r*rsc;

    Cc=1./(1+Rc.*Ra./(Rs.*(Rc+Ra)));
    Cs=1./(1+Rs.*Ra./(Rc.*(Rs+Ra)));

    [fenzi1]=slp.*(Rn-G)+(d_a*Cp*VPD-slp.*rac.*(Rns-G))./(raa+rac);
    [fenmu1]=slp+r*(1+rsc./(raa+rac));

    [PMc]=fenzi1./fenmu1;

    [fenzi2]=slp.*(Rn-G)+(d_a*Cp*VPD-slp.*ras.*(Rn-Rns))./(raa+ras);
    [fenmu2]=slp+r*(1+rss./(raa+ras));
    [PMs]=fenzi2./fenmu2;

    [T]=Cc.*PMc/2450; % 1wa/s=1 J, 1 g H2O=2450 J
    [E]=Cs.*PMs/2450;
    [ET]=T+E;
    [E_ET]=E./ET;

    shuttle(:,21)=T;% from g s-1 to kg day-1
    shuttle(:,22)=E;
    shuttle(:,23)=ET;
    shuttle(:,24)=E_ET;
    shuttle(:,25)=rsc;
    shuttle(:,26)=rss;
    shuttle(:,27)=raa;
    % decoupling coeficient
   % shuttle(:,26)=(slp/r+1)./(slp/r+1+(rsc./(raa+ras)));




