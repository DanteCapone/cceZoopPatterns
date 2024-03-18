addpath 'm_map1.4'

close all
clear all

%Choose boundaries
%latmin=32;
latmin=33;% this is for AT
%latmax=38;
%latmax=36; %this is for AT
latmax=37;
lonmin=-116;
%lonmax=-123; %this is for AT
%lonmax=-129; %this is for CCT
lonmax=-131.5; %this is to include the offshore

%Reading in coastline, CalCOFI station, and bathymetry data
load('Data\NEPac_MedRes_Coastline.mat')

CalCOFI=load('Data\CCStaPos.csv');

load('Data\Global4min.mat')
lat=-lat;
bathymetry=z;
bathymetry(find(bathymetry>0))=4000;
bathymetry=bathymetry';



z(find(z>0))=z(find(z>0))/2+2500;

%CCE-P1908 drifter data
P1908array=xlsread('Data\CCE-P1908 - DrifterTracks.xls','DrifterTracks','H2:I1870');
temp=xlsread('Data\CCE-P1908 - DrifterTracks.xls','DrifterTracks','C2:C1870');
P1908array(:,3)=temp;
P1908sedtrap=xlsread('Data\CCE-P1908 - DrifterTracks.xls','DrifterTracks','H1871:I3887');
temp=xlsread('Data\CCE-P1908 - DrifterTracks.xls','DrifterTracks','C1871:C3887');
P1908sedtrap(:,3)=temp;
ZoopB=xlsread('bridgelog_draftReport.xlsx','k2:l50');
ZoopMoc=xlsread('bridgelog_draftReport.xlsx','ag2:ah30');
ZoopRing=xlsread('bridgelog_draftReport.xlsx','bb2:bc50');
ZoopSalp=xlsread('bridgelog_draftReport.xlsx','bx2:by50');

%% Reading in Met Data
[shippos,temperature,salinity,fluorescence,flowmeter] = GetMetData(0715,0811);
salinity(find(salinity<30))=NaN;
fluorescence(find(fluorescence<0.05))=NaN;



%Reading in realtime drifter data
startDate = '2021-06-18';
startTime = '0130:00';endTime = datestr(now,13);
endTime = [endTime(1:2) endTime(4:8)];
endDate = datestr(now,29);
[deviceNumPG,serialTimePG,datLatPG,datLonPG,battVoltPG,gpsQualPG,subFlagPG,...
        pos1PG,pos3PG,pos4PG,time1PG,time3PG,time4PG] = getLatestWebData(startDate,startTime,endDate,endTime);
SedTrap = [datLatPG(find(deviceNumPG == 3)),datLonPG(find(deviceNumPG == 3)),serialTimePG(find(deviceNumPG == 3))];
Array = [datLatPG(find(deviceNumPG == 4)),datLonPG(find(deviceNumPG == 4))];
%SedTrap_degmin = [floor(SedTrap(:,1)),(SedTrap(:,1)-floor(SedTrap(:,1)))*60



%----------------------------Plotting Fig 1 (SST)--------------
figure('Position',[100 -100 1400 1200])
hold on
m_proj('lambert','lon',[lonmin lonmax],'lat',[latmin latmax]);

m_scatter(shippos(:,1),shippos(:,2),60,temperature);
%m_pcolor(lon,lat,double(bathymetry))
%shading interp
m_plot(CalCOFI(:,2),CalCOFI(:,3),'oy','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',4)

%m_plot(P1908sedtrap(:,2),P1908sedtrap(:,1),'.','Color','r')
%m_plot(P1908array(:,2),P1908array(:,1),'.','Color','y')
%m_plot(SedTrap(:,2),SedTrap(:,1),'or','MarkerFaceColor','r')
%m_plot(Array(:,2),Array(:,1),'m','LineWidth',3)
%m_plot(SedTrap(:,2),SedTrap(:,1),'*k','MarkerFaceColor','k')
%m_plot(Array(:,2),Array(:,1),'k','LineWidth',3)
Mkp=10;l=1.5;
b = m_plot(ZoopB(2:end,2),ZoopB(2:end,1),'ok','MarkerFaceColor','w', 'MarkerEdgeColor', 'k',...
    'Markersize', Mkp, 'linewidth', l);
m = m_plot(ZoopMoc(2:end,2),ZoopMoc(2:end,1),'sk','MarkerFaceColor','w', 'MarkerEdgeColor', 'k',...
    'Markersize', Mkp, 'linewidth', l);
%r = m_plot(ZoopRing(:,2),ZoopRing(:,1),'*k','MarkerFaceColor','k');
%s = m_plot(ZoopSalp(:,2),ZoopSalp(:,1),'hk','MarkerFaceColor','w');
%b1 = m_plot(ZoopB(32:39,2),ZoopB(32:39,1),'ok','MarkerFaceColor','k', 'Markersize', 6);
%b2 = m_plot(ZoopB(40:47,2),ZoopB(40:47,1),'ok','MarkerFaceColor','w', 'Markersize', 6, ...
%    'Linewidth', 2);
m_plot(ZoopB(43,2),ZoopB(43,1),'ok','MarkerFaceColor','k', 'Markersize', 2);

F=15;mk=8;
m_plot(NEPac_MedRes_Coastline(:,1),NEPac_MedRes_Coastline(:,2),'k','LineWidth',4)
%m_plot(-117.1611,32.7157,'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
%m_text(-117.1611+0.1,32.7157+0.3,['San',char(10),'Diego'],'FontSize',F)
m_plot(-120.4716,34.4486,'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',mk)
m_text(-120.4716+0.05,34.4486+0.2,['Point Conception'],'FontSize',F)
%m_plot(-121.8947,36.6002,'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',mk)
%m_text(-121.8947+0.1,36.6002,['Monterey'],'FontSize',F)
m_plot(-120.8500,35.3659,'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',mk)
m_text(-120.8500+0.05,35.3659+0.1,['Morro Bay'],'FontSize',F)
m_plot(-121.9016,36.3064,'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',mk)
m_text(-121.9016+0.05,36.3064+0.05,['Point Sur'],'FontSize',F)
m_grid
%legend([ b m r s], {'Bongo', 'MOCNESS', 'Ring net', 'Salp net'})
legend([ b m], {'Bongo', 'MOCNESS'})
%legend([ b1 b2], {'California Current Transect (CCT)', 'Alongshore Transect (AT)'}, 'Fontsize', F, 'Location', 'NW')
hold on

%colormap(winter)
colormap (jet)
%colorbar
h = colorbar;
set(get(h,'label'),'string','Temperature (°C)', 'Fontsize', F);
%title('SST')

