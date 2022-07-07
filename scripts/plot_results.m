clear;

%% Load run data

% Meltwater Output Directory
outDirectory='processed/';

% Specify Run to Process
runDate='20211115_ADJ_M3/';
runname= 'basin-multi-adj-ekh-alb-2007-streams.mat';

% Load data
path2output=[outDirectory runDate runname];
load(path2output);

years2run = (1996:2013);
nyears = length(years2run);

%% Re-order total runoff volume matrix

% basinOrder = [10 21 29 34 41 43 45 61 62 64 65 66 50 71 74];

% Basin ID
basinkey=[...
    62,41,45,74,66,...
    50,43,65,34,61,...
    29,71,21,10,64];

% Remove santa fe
basinOrder = [21 29 34 41 43 45 61 62 64 65 66 50 71 74];

for b=1:14
    doB = find(basinkey == basinOrder(b));
    for yr=years2run
        streamYrVoltotal(yr-1995,b)=streamYrVol(yr-1995,doB);
        modelYrVoltotal(yr-1995,b)=modelSmVol(yr-1995,doB);
    end
end

for b=1:14
    doB = find(basinkey == basinOrder(b));
    for yr=years2run
        streamYrVolpercent(yr-1995,b)=streamYrVol(yr-1995,doB)/sum(streamYrVoltotal(yr-1995,:))*100;
        modelYrVolpercent(yr-1995,b)=modelSmVol(yr-1995,doB)/sum(modelYrVoltotal(yr-1995,:))*100;
    end
end

streamYrLaketotal = zeros(18,4);
modelSmLaketotal = zeros(18,4);
modelMeasLaketotal = zeros(18,4);

for b=1:14
    doB = find(basinkey == basinOrder(b));
    if b < 3
        lake = 1;
    elseif b == 3 | b == 4
        lake = 2;
    elseif b > 4 & b < 14
        lake = 3;
    else 
        lake = 4;
    end
    for yr=years2run
        streamYrLaketotal(yr-1995,lake)= streamYrLaketotal(yr-1995,lake) + streamYrVol(yr-1995,doB);
        modelSmLaketotal(yr-1995,lake)= modelSmLaketotal(yr-1995,lake) + modelSmVol(yr-1995,doB);
        modelMeasLaketotal(yr-1995,lake)= modelMeasLaketotal(yr-1995,lake) + modelYrVolMeas(yr-1995,doB);
    end
end

for lake=1:4
    for yr=years2run
        streamYrLakepercent(yr-1995,lake)=streamYrLaketotal(yr-1995,lake)/sum(streamYrLaketotal(yr-1995,:))*100;
        modelSmLakepercent(yr-1995,lake)=modelSmLaketotal(yr-1995,lake)/sum(modelSmLaketotal(yr-1995,:))*100;
        modelMeasLakepercent(yr-1995,lake)=modelMeasLaketotal(yr-1995,lake)/sum(modelMeasLaketotal(yr-1995,:))*100;
    end
end

%% Lake Basins - Total runoff volume per season by lake - LINE PLOT

% Create plot
figure(10001); clf; clear ha; ha = tight_subplot(3,1, [0.03 0.05], [.11 .01], [.08 .03]);
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.8 0.9])

x=[years2run(2):years2run(end)];

myC= [37 52 148 % Bonney
    239 59 44 % Canada
    0 69 41 % Fryxell
    70 75 85]; % McMurdo

axes(ha(1));

% Measured streamflow Bonney
Ymeas=sum(streamYrLaketotal(2:nyears,1),2)/1000000;
err=Ymeas*0.2;
errorbar(x,Ymeas,err, '-k', 'LineWidth', 1.5); hold on;

% M4 Measured Bonney
% plot(x,(sum(modelMeasLaketotal(1:nyears,1),2)/1000000)', ':','Color', (myC(1,:)/255), 'LineWidth', 2);
plot(x,(sum(modelMeasLaketotal(2:nyears,1),2)/1000000)', '-.k','LineWidth', 1.5);
% plot(x,(sum(modelSmLaketotal(1:nyears,1),2)/1000000)', '-.k','Color', (myC(1,:)/255), 'LineWidth', 1);

% ylim([0 5])
% xlim([1995 2014])
% set(gca,'xtick',(1995:1:2014));
xlim([1996 2014])
set(gca,'xtick',(1996:1:2014));
set(gca,'xticklabel',' ');
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize',14,'LineWidth', 1.5, 'GridColor', 'k');

axes(ha(2));

% Measured streamflow Hoare
Ymeas=sum(streamYrLaketotal(2:nyears,2),2)/1000000;
err=Ymeas*0.2;
errorbar(x,Ymeas,err, '-k', 'LineWidth', 1.5); hold on;

% M4 Measured Hoare
% plot(x,(sum(modelMeasLaketotal(1:nyears,2),2)/1000000)', ':','Color', (myC(2,:)/255), 'LineWidth', 2);
plot(x,(sum(modelMeasLaketotal(2:nyears,2),2)/1000000)', '-.k','LineWidth', 1.5);
% plot(x,(sum(modelSmLaketotal(1:nyears,2),2)/1000000)', '--','Color', (myC(2,:)/255), 'LineWidth', 1);

% ylim([0 5])
% xlim([1995 2014])
% set(gca,'xtick',(1995:1:2014));
xlim([1996 2014])
set(gca,'xtick',(1996:1:2014));
set(gca,'xticklabel',' ');
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize',14,'LineWidth', 1.5, 'GridColor', 'k');

axes(ha(3));

% Measured streamflow Fryxell
Ymeas=sum(streamYrLaketotal(2:nyears,3),2)/1000000;
err=Ymeas*0.2;
errorbar(x,Ymeas,err, '-k', 'LineWidth', 1.5); hold on;

% M4 Measured Fryxell
% plot(x,(sum(modelMeasLaketotal(1:nyears,3),2)/1000000)', ':','Color', (myC(3,:)/255), 'LineWidth', 2);
plot(x,(sum(modelMeasLaketotal(2:nyears,3),2)/1000000)', '-.k', 'LineWidth', 1.5);
% plot(x,(sum(modelSmLaketotal(1:nyears,3),2)/1000000)', '--','Color', (myC(3,:)/255), 'LineWidth', 1);

% ylim([0 5])
% xlim([1995 2014])
% set(gca,'xtick',(1995:1:2014));
xlim([1996 2014])
set(gca,'xtick',(1996:1:2014));
xtickangle(45)
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize',14,'LineWidth', 1.5, 'GridColor', 'k');

axes(ha(2));
ylabel('Annual Discharge [x10^{6} m^3 yr^{-1}]', 'FontSize', 14)
set(gca,'yTickLabelMode','auto')

axes(ha(3));
xlabel('Water Year', 'FontSize', 14)

% Legend
legendArray = {['Obs. ' char(177) '20%'],...
    'Sim. M4 - days with obs.', 'Sim. M4 - all days', 'Sim. M4 - days with obs.'};

legend(legendArray, 'Location', 'northwest');
legend boxoff

%% Model vs. measured runoff daily - 1:1 SCATTERPLOTS All

fig = figure (100); clf; hold on;
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.5 0.75]);

% Loop through basins all years
for b=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]
    
    % Plot 1:1 Scatter Plot
    x = streams.TDQ(:,b);
    y = modelTDQ(:,b);
    sz = 10;
    c = day(datetime(streams.dates(:,:), 'convertfrom', 'datenum')+caldays(180),'dayofyear');
    
    scatter(x,y,sz,c,'filled')
    colorbar()
    caxis([150, 200])
    
    line([1 100000], [1 100000],'Color','red','LineStyle','--','LineWidth', 1.5)
    
    set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize',14,'LineWidth', 1.5, 'GridColor', 'k', 'box', 'on');
    set(gca,'XScale','log','XLim',[1 100000])
    set(gca,'YScale','log','YLim',[1 100000])
    xlabel('Log of Obs. Discharge [m^3 d^{-1}]', 'FontSize', 18)
    ylabel('Log of Sim. Discharge [m^3 d^{-1}]', 'FontSize', 18)
end

%% Model vs. measured runoff annual - 1:1 SCATTERPLOTS All

fig = figure (101); clf; hold on;
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.5 0.75]);

cmap = parula(15);
markers = ['+','o','*','x','v','d','^','s','>','<','+','o','*','x','v'];

% Loop through basins all years
for b=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]
    
    % Plot 1:1 Scatter Plot
    x = streamYrVol(:,b);
    y = modelYrVolMeas(:,b);
    sz = 50;
    c = cmap(b,:);
    mkr = markers(b);
    
    
    scatter(x,y,sz,c,mkr)
    
%     colorbar('YTickLabel',basinkey);
    caxis([0 16]);
    
    line([1 10000000], [1 10000000],'Color','red','LineStyle','--','LineWidth', 1.5, 'HandleVisibility','off')
    
    set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize',14,'LineWidth', 1.5, 'GridColor', 'k', 'box', 'on');
    set(gca,'XScale','log','XLim',[100 10000000])
    set(gca,'YScale','log','YLim',[100 10000000])
    xlabel('Log of Obs. Discharge [m^3 d^{-1}]', 'FontSize', 18)
    ylabel('Log of Sim. Discharge [m^3 d^{-1}]', 'FontSize', 18)
end

legend(streamname)
%% Lake Basins - Percent contribution to total discharge - BAR CHART


myC= [37 52 148 % Bonney
    239 59 44 % Canada
    0 69 41 % Fryxell
    70 75 85]; % McMurdo

figure(10004); clf; clear ha;  ha = tight_subplot(2,2, [0.1 0.07], [.1 .04], [.1 .06]);
set(gcf, 'Name', 'all-season-contribution');
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.8])

% Measured discharge absolute
% subplot(2,2,1)
axes(ha(1));
bar1 = bar(years2run, streamYrLaketotal, 1,'stack');
for l=1:4
    set(bar1(l),'facecolor',(myC(l,:)/255))
end
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'GridColor', 'k');

ylabel('Annual Discharge [m^3 yr^{-1}]', 'FontSize', 14)
set(gca,'yTickLabelMode','auto')
axis([1994 2013 0 8000000]);

% Measured discharge percent
% subplot(2,2,2)
axes(ha(3));
bar2 = bar(years2run, streamYrLakepercent, 1,'stack');
for l=1:4
    set(bar2(l),'facecolor',(myC(l,:)/255))
end
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'GridColor', 'k');

ylabel('% of Total Annual Discharge', 'FontSize', 14)
set(gca,'yTickLabelMode','auto')
set(gca,'yticklabel',num2str(get(gca,'ytick')'))  % convert to 1000 m^3 units

xlabel('Season', 'FontSize', 14)
axis([1994 2013 0 100]);

% Modeled discharge absolute
% subplot(2,2,3)
axes(ha(2));
bar3 = bar(years2run, modelMeasLaketotal, 1,'stack');
for l=1:4
    set(bar3(l),'facecolor',(myC(l,:)/255))
end
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'GridColor', 'k');

%ylabel('Annual Discharge [m^3 yr^{-1}]', 'FontSize', 12)
set(gca,'yTickLabelMode','auto')
axis([1994 2013 0 8000000]);

legendArray = {'Lake Bonney', 'Lake Hoare', 'Lake Fryxell', 'McMurdo Sound'};
legend(legendArray,'Location','Best');
legend boxoff

% Modeled discharge percent
% subplot(2,2,4)
axes(ha(4));
bar4 = bar(years2run, modelMeasLakepercent, 1,'stack');
for l=1:4
    set(bar4(l),'facecolor',(myC(l,:)/255))
end
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'GridColor', 'k');

set(gca,'yTickLabelMode','auto')
xlabel('Season', 'FontSize', 14)
axis([1994 2013 0 100]);

%% Streams - Percent contribution to total discharge - BAR CHART

figure(10005); clf; clear ha;  ha = tight_subplot(2,2, [0.1 0.07], [.1 .04], [.1 .06]);
set(gcf, 'Name', 'all-season-contribution');
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.8])

% Colors and Legend
for b=1:14
    streamnameOrder(b,1) = streamname(find(basinkey == basinOrder(b)));
end

myC= [34 94 168 %Priscu
    29 145 192 %Lawson
    239 59 44 %House
    251 106 74 %Andersen
    0 69 41 %Green
    0 104 55 %Canada
    35 132 67 %Delta
    65 171 93 %Crescent
    120 198 121 %Huey
    173 221 142 %Aiken
    217 240 163 %VonGuerard
    247 252 185 %Harnish
    255 255 229 %Lost Seal
    70 75 85]; %Commonwealth

% myC= [37 52 148 %Santa Fe
%     34 94 168 %Priscu
%     29 145 192 %Lawson
%     239 59 44 %House
%     251 106 74 %Andersen
%     0 69 41 %Green
%     0 104 55 %Canada
%     35 132 67 %Delta
%     65 171 93 %Crescent
%     120 198 121 %Huey
%     173 221 142 %Aiken
%     217 240 163 %VonGuerard
%     247 252 185 %Harnish
%     255 255 229 %Lost Seal
%     70 75 85]; %Commonwealth

% Measured discharge absolute
% subplot(2,2,1)
axes(ha(1));
bar1 = bar(years2run, streamYrVoltotal, 1,'stack');
for k=1:14
    set(bar1(k),'facecolor',(myC(k,:)/255))
end
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'GridColor', 'k');

ylabel('Annual Discharge [m^3 yr^{-1}]', 'FontSize', 14)
set(gca,'yTickLabelMode','auto')
axis([1994 2013 0 8000000]);

% Measured discharge percent
% subplot(2,2,2)
axes(ha(3));
bar2 = bar(years2run, streamYrVolpercent, 1,'stack');
for k=1:14
    set(bar2(k),'facecolor',(myC(k,:)/255))
end
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'GridColor', 'k');

ylabel('% of Total Annual Discharge', 'FontSize', 14)
set(gca,'yTickLabelMode','auto')
set(gca,'yticklabel',num2str(get(gca,'ytick')'))  % convert to 1000 m^3 units

xlabel('Season', 'FontSize', 14)
axis([1994 2013 0 100]);

% Modeled discharge absolute
% subplot(2,2,3)
axes(ha(2));
bar3 = bar(years2run, modelYrVoltotal, 1,'stack');
for k=1:14
    set(bar3(k),'facecolor',(myC(k,:)/255))
end
axis([1994 2013 0 8]);
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'GridColor', 'k');

%ylabel('Annual Discharge [m^3 yr^{-1}]', 'FontSize', 12)
set(gca,'yTickLabelMode','auto')
axis([1994 2013 0 8000000]);

% legendOrder = [15 14 13 12 11 10 9 8 7 6 5 4 3 2 1];
legendOrder = [14 13 12 11 10 9 8 7 6 5 4 3 2 1];
legend(fliplr(bar3),streamnameOrder(legendOrder),'Location','Best', 'FontSize', 8);
legend boxoff

% Modeled discharge percent
% subplot(2,2,4)
axes(ha(4));
bar4 = bar(years2run, modelYrVolpercent, 1,'stack');
for k=1:14
    set(bar4(k),'facecolor',(myC(k,:)/255))
end

set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'GridColor', 'k');

set(gca,'yTickLabelMode','auto')
xlabel('Season', 'FontSize', 14)
axis([1994 2013 0 100]);

%% Taylor Valley - Model vs. measured runoff - SEASONAL HYDROGRAPH

figure(1111); clf;
set(gcf,'units','normalized','outerposition',[0.1 0.3 0.8 0.4])

yr = 2001;
b = 9;

strt=find(streams.dates==datenum([yr 11 1]));
fnsh=find(streams.dates==datenum([yr+1 3 1]));

confplot(streams.dates(strt:fnsh),streams.TDQ(strt:fnsh,b),...
    streams.TDQ(strt:fnsh,b)*0.20,streams.TDQ(strt:fnsh,b)*0.20,'-k','LineWidth', 1);
hold on
plot(streams.dates(strt:fnsh),modelTDQ(strt:fnsh,b), '-.k', 'LineWidth', 1.5);
% plot(streams.dates(strt:fnsh),run2.modelTDQ(strt:fnsh,b), ':k', 'LineWidth', 2);

xlim([streams.dates(strt) streams.dates(fnsh)]);

datetick('x','keeplimits');
%ylim([0 max([streams.TDQ(:,b)])*2]);
ylabel('Discharge [m^3 d^{-1}]');
yticks('auto')
xlabel('Date', 'FontSize', 14)

%title([streamname(b), num2str(yr), 'to', num2str(yr+1)]);
set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize',14,'LineWidth', 1.5, 'GridColor', 'k');

%% Cummulative Error

% load('modelSkill_Tmp.mat');
% 
% cumError = cumsum(err);
% 
% % Create plot
% figure (6); clf; hold all;
% set(gcf, 'Name', 'all-volume-compare');
% set(gcf,'units','normalized','outerposition',[0.2 0.4 0.6 0.5])
% set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'GridColor', 'k');
% 
% % Absolute
% yyaxis left
% plot(streams.dates,cumError/1000000,'-k','LineWidth', 1);
% ylim([-1 25])
% ylabel('Error (10^6 m^3)', 'FontSize', 18)
% 
% xlabel('Year', 'FontSize', 18)
% xlim([streams.dates(1)-(2*365) streams.dates(end)]);
% datetick('x','keeplimits');
% xtickangle(45)
% 
% % Normalized
% yyaxis right
% plot(streams.dates,cumError/13650000,'--k','LineWidth', 1);
% ylim([-0.075 1.8315])
% ylabel('Area Normalized Error (m)', 'FontSize', 18)
% 
% xlabel('Year', 'FontSize', 18)
% xlim([streams.dates(1)-(2*365) streams.dates(end)+(2*365)]);
% datetick('x','keeplimits');
% xtickangle(45)
% 
% legendArray = {'Cummulative Error (Obs-Mod)'};
% %legendArray = {'Q_{OBS}','Q_{MOD}', 'Cummulative Error (Obs-Mod)'};
% legend(legendArray, 'Location', 'southeast');
% legend boxoff

%% Taylor Valley - Model vs. measured runoff - HYDROGRAPHS

% years2run = [1995 2001 2008 2010];
% 
% % Loop through basins
% for b=[1 2 3 4 5 6 7 8 9 10]
%     fig = figure (b+100); clf;
%     set(gcf, 'Name', char(streamname(b)));
%     set(gcf,'units','normalized','outerposition',[0 0.0375 0.5 0.96]);
% 
%     % Loop through years
%     i = 1;
%     for yr=years2run
%         
%         nrows = 4;
%         ncols = 1;
%         
%         % Plot subplot
%         subplot(nrows,ncols,i);
%         strt=find(streams.dates==datenum([yr 11 1]));
%         fnsh=find(streams.dates==datenum([yr+1 3 1]));
%         
%         confplot(streams.dates(strt:fnsh),streams.TDQ(strt:fnsh,b),...
%             streams.TDQ(strt:fnsh,b)*0.20,streams.TDQ(strt:fnsh,b)*0.20,'-k','LineWidth', 1);
%         hold on
% %         plot(streams.dates(strt:fnsh),streams.ErrUpr(strt:fnsh,b), '-r', 'LineWidth', 1);
% %         plot(streams.dates(strt:fnsh),streams.ErrLwr(strt:fnsh,b), '-b', 'LineWidth', 1);
%         
%         plot(streams.dates(strt:fnsh),modelTDQ(strt:fnsh,b), '--k', 'LineWidth', 1.5);
%         
%         xlim([streams.dates(strt) streams.dates(fnsh)]);
%         
%         datetick('x','keeplimits');
%         %ylim([0 max([streams.TDQ(:,b)])*2]);
%         ylabel('Discharge [m^3 d^{-1}]');
%         yticks('auto')
%         xlabel('Date', 'FontSize', 14)
%         set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize',14,'LineWidth', 1.5, 'GridColor', 'k');
%         
%         i = i + 1;
%     end
% end

%% Model vs. measured runoff - 1:1 SCATTERPLOTS All One Plot by Basin

fig = figure (100); clf; clear ha;  ha = tight_subplot(5,3, [.04 .04], [.07 .04], [.1 .2]);
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);

i = 1;

basinOrder = [10 21 29 34 41 43 45 50 66 61 62 64 65 71 74];

% Loop through basins all years
for b=1:15
    
    doB = find(basinkey == basinOrder(b));
    
    axes(ha(i)); hold all;
    title(char(streamname(doB)));
    
    for yr=1995:2011
        
        strt=find(streams.dates==datenum([yr 11 1]));
        fnsh=find(streams.dates==datenum([yr+1 3 1]));
        
        % Plot 1:1 Scatter Plot
        x = streams.TDQ(strt:fnsh,doB);
        y = modelTDQ(strt:fnsh,doB);
        sz = 25;
        c = datenum([1995 11 1]) + streams.dates(strt:fnsh,1) - datenum([yr 11 1]);
        
        scatter(x,y,sz,c,'filled'); hold on;
        line([1 100000], [1 100000],'Color','red','LineStyle','--','LineWidth', 1.5)
        
    end
    
    set(gca,'XColor','k', 'YColor', 'k','LineWidth', 1.5, 'FontWeight', 'bold', 'GridColor', 'k', 'box', 'on')
    set(gca,'XScale','log','XLim',[1 100000])
    set(gca,'YScale','log','YLim',[1 100000])
    set(gca,'xticklabelmode','auto')
    set(gca,'yticklabelmode','auto')
   
    i = i + 1;
    
end

axes(ha(7))
ylabel('Log of Sim. Discharge [m^3 d^{-1}]', 'FontSize', 18)

axes(ha(14))
xlabel('Log of Obs. Discharge [m^3 d^{-1}]', 'FontSize', 18)

axes(ha(15))
colorbar()
set(gca, 'FontSize', 14)
hp = get(gca,'Position');
colorbar('Position', [hp(1)+hp(3)+0.05  hp(2)  0.05  hp(2)+hp(3)*4])
cbdate({'Nov 1, 95';'Dec 1, 95';'Jan 1, 96';'Feb 1, 96';'Mar 1, 96'})
cbdate('dd-mmm.')

%% Model vs. measured runoff - 1:1 SCATTERPLOTS All One Plot by Year

fig = figure (100); clf; clear ha;  ha = tight_subplot(5,3, [.04 .04], [.07 .04], [.1 .2]);
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);

colormap(winter)

i = 1;

basinOrder = [10 21 29 34 41 43 45 50 66 61 62 64 65 71 74];

% Loop through years all basins
for b=1:15
    
    doB = find(basinkey == basinOrder(b));
    
    axes(ha(i)); hold all;
    title(char(streamname(doB)));
    
    % Plot 1:1 Scatter Plot
    x = streams.TDQ(:,doB);
    y = modelTDQ(:,doB);
    sz = 25;
    c = streams.dates(:,1);
    
    scatter(x,y,sz,c,'filled'); hold on;
    line([1 100000], [1 100000],'Color','red','LineStyle','--','LineWidth', 1.5)
    
    
    set(gca,'XColor','k', 'YColor', 'k','LineWidth', 1.5, 'FontWeight', 'bold', 'GridColor', 'k', 'box', 'on')
    set(gca,'XScale','log','XLim',[1 100000])
    set(gca,'YScale','log','YLim',[1 100000])
    set(gca,'xticklabelmode','auto')
    set(gca,'yticklabelmode','auto')
    
    i = i + 1;
    
end

axes(ha(7))
ylabel('Log of Sim. Discharge [m^3 d^{-1}]', 'FontSize', 18)

axes(ha(14))
xlabel('Log of Obs. Discharge [m^3 d^{-1}]', 'FontSize', 18)

axes(ha(15))
colorbar()
hp = get(gca,'Position');
colorbar('Position', [hp(1)+hp(3)+0.02  hp(2)  0.05  hp(2)+hp(3)*4], 'FontSize', 14)
cbdate('mmm. yyyy')

%% Model vs. measured runoff - 1:1 SCATTERPLOTS All One Plot by Year by Basin

% % Loop through basins loop through years
% for yr=1995:2011
%     
%     fig = figure (yr); clf;
%     set(gcf,'units','normalized','outerposition',[0.01 0.01 0.9 0.9]);
%     
%     i = 1;
%     
%     strt=find(streams.dates==datenum([yr 11 1]));
%     fnsh=find(streams.dates==datenum([yr+1 3 1]));
%     
%     for b=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]
%         
%         subplot(3,5,i);
%         
%         % Plot 1:1 Scatter Plot
%         x = streams.TDQ(strt:fnsh,b);
%         y = modelTDQ(strt:fnsh,b);
%         sz = 25;
%         c = streams.dates(strt:fnsh,1);
%         
%         scatter(x,y,sz,c,'filled')
%         hold all;
%         set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize',14,'LineWidth', 1.5, 'GridColor', 'k', 'box', 'on');
%         set(gca,'XScale','log','XLim',[1 100000])
%         set(gca,'YScale','log','YLim',[1 100000])
%         
%         i = i +1;
%         
%     end
%     
%     colorbar()
%     cbdate()
%     
% end
%% Model vs. measured runoff - Hydrographs and 1:1 SCATTERPLOTS all years

% yr = 1995;
% 
% % Loop through basins
% for b=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]
%     fig = figure (b+200); clf;
%     set(gcf, 'Name', char(streamname(b)));
%     set(gcf,'units','normalized','outerposition',[0 0.5 0.75 0.5]);
%     
%     i = 1;
%     
%     nrows = 1;
%     ncols = 3;
%     
%     % Plot subplot
%     subplot(nrows,ncols,i:i+1);
%     strt=find(streams.dates==datenum([yr 11 1]));
%     fnsh=find(streams.dates==datenum([yr+17 3 1]));
%     
%     confplot(streams.dates(strt:fnsh),streams.TDQ(strt:fnsh,b),...
%         streams.TDQ(strt:fnsh,b)*0.20,streams.TDQ(strt:fnsh,b)*0.20,'-k','LineWidth', 1);
%     hold on
%     plot(streams.dates(strt:fnsh),streams.ErrUpr(strt:fnsh,b), '-r', 'LineWidth', 1);
%     plot(streams.dates(strt:fnsh),streams.ErrLwr(strt:fnsh,b), '-b', 'LineWidth', 1);
%     
%     plot(streams.dates(strt:fnsh),modelTDQ(strt:fnsh,b), '--k', 'LineWidth', 2);
% %     plot(streams.dates(strt:fnsh),run2.modelTDQ(strt:fnsh,b), ':k', 'LineWidth', 2);
%     
%     xlim([streams.dates(strt) streams.dates(fnsh)]);
%     
%     datetick('x','keeplimits');
%     %ylim([0 max([streams.TDQ(:,b)])*2]);
%     ylabel('Discharge [m^3 d^{-1}]');
%     yticks('auto')
%     xlabel('Date', 'FontSize', 14)
%     set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize',14,'LineWidth', 1.5, 'GridColor', 'k');
%     
%     % Plot 1:1 Scatter Plot
%     subplot(nrows,ncols,i+2);
%     scatter(streams.TDQ(strt:fnsh,b),modelTDQ(strt:fnsh,b), 10,'k','filled');
%     
%     %title([streamname(b), num2str(yr), 'to', num2str(yr+1)]);
%     set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize',14,'LineWidth', 1.5, 'GridColor', 'k', 'box', 'on');
%     
%     i = i + 3;
% end

%% Taylor Valley - Total runoff volume per season by lake - LINE and BAR PLOT

% % Create plot
% figure (1); clf; hold all;
% set(gcf, 'Name', 'all-volume-compare');
% set(gcf,'units','normalized','outerposition',[0.2 0.4 0.6 0.5])
% set(gca,'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'FontSize', 10, 'GridColor', 'k');
% 
% % Measured lake inflow
% hold on
% %bar(years2run, streamYrLaketotal/1000000, 1,'stack');
% % Modeled lake inflow
% bar(years2run, modelSmLaketotal/1000000, 1,'stack');
% 
% % measured runoff
% plot(years2run,sum(streamYrVol(1:nyears,:)/1000000,2),'-k','LineWidth', 1);
% ylim([0 8000000/1000000])
% % modeled runoff
% plot(years2run, sum(modelSmVol(1:nyears,:),2)/1000000, 's--k', 'LineWidth', 1, 'MarkerSize', 6)
% 
% ylabel('Glacier Runoff (10^6 m^6 yr^{-1})', 'FontSize', 9)
% set(gca,'yTickLabelMode','auto')
% set(gca,'yticklabel',num2str(get(gca,'ytick')'))  % convert to 1000 m^3 units
% 
% xlabel('Melt Season', 'FontSize', 9)
% xlim([1994 2013])
% set(gca,'xtick',(1994:1:2013));
% xtickangle(45)
% 
% legendArray = {'Lake Bonney', 'Lake Hoare', 'Lake Fryxell', 'McMurdo Sound', 'Q_{OBS}','Q_{MOD1}'};
% legend(legendArray, 'Location', 'best');
% legend boxoff
