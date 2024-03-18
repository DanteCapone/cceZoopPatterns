%% CCE Physical Maps

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


figure('Position',[100 -100 1400 1200])
hold on
m_proj('lambert','lon',[lonmin lonmax],'lat',[latmin latmax]);

m_scatter(met.shippos(:,1),met.shippos(:,2),60,met.temperature);
m_coast('patch',[.9 .9 .9],'edgecolor','none');
%%

axes('position',[.55 .35 .37 .37]);
    m_proj('albers equal-area','lat',[latmin latmax],'long',[lonmin lonmax],'rect','on');
    m_gshhs_f('patch',[.7 .9 .7]);
    m_grid('linestyle','none','linewidth',2,'tickdir','out',...
           'xaxisloc','top','yaxisloc','right','fontsize',6);
    m_text(-63.95,46.56,'GSHHS\_F (full)','color','m','fontweight','bold');
    m_ruler([.5 .8],.2,3,'fontsize',8);