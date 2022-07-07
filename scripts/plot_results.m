clear;

%% Load run data

% Meltwater Output Directory
outDirectory='../processed/';

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

%% Single Stream - Model vs. measured runoff - SEASONAL HYDROGRAPH

figure(2); clf;
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

%% Lake Basins - Total runoff volume per season by lake - LINE PLOT

% Create plot
figure(3); clf; clear ha; ha = tight_subplot(3,1, [0.03 0.05], [.11 .01], [.08 .03]);
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