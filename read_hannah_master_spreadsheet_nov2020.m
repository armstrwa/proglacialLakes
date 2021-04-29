%% Script to read Hannah's master lake spreadsheet
%  Also generally useful for reading in csv's with multiple data types
clear all
close all
clc

%%
% data file
%filename = '/Users/williamarmstrong/Google Drive/appstate/projects/lakes/data/Master_Lake_Data_20200205.csv';
%filename = '/Users/williamarmstrong/Google Drive/appstate/projects/lakes/data/Master_Lake_Data_Numerical_200717.csv';
%path = 'C:\Users\armstrongwh\Google Drive\appstate\projects\lakes\data\';%lab cpu
%path = '/Users/williamarmstrong/Google Drive/appstate/projects/lakes/manuscripts/papers/hannahPaperFinal/';
path = 'C:\Users\armstrongwh\Google Drive\appstate\projects\lakes\manuscripts\papers\hannahPaperFinal\'
%path = '/Users/williamarmstrong/Google Drive/appstate/projects/lakes/data/'; % home cpu
fn = 'Master_Lake_Data_Numerical_corrected_whaUpdateWRelAreaAndCoastDist.csv';
filename = [path fn];

% - Get structure from first line.
fid  = fopen( filename, 'r' ) ;
% do this to skip multiple header lines
line = fgetl( fid ) ;
line = fgetl( fid ) ;
line = fgetl( fid ) ;
fclose( fid ) ;
%isStrCol = isnan( str2double( regexp( line, '[^\t]+', 'match' ))) ;
isStrCol = isnan( str2double( regexp( line, '[^,]+', 'match' ))) ; % this assumes a comma delimiter
%isStrCol(108) = 1; % col 108 (second associated glacier ID) was breaking because first were -999 (if none exist), and then at 34th lake, has string (the associated RGI ID)
% - Build formatSpec for TEXTSCAN.
fmt = cell( 1, numel(isStrCol) ) ;
fmt(isStrCol)  = {'%s'} ; % string format
fmt(~isStrCol) = {'%f'} ; % floating point number format
fmt = [fmt{:}] ;
% - Read full file.
fid  = fopen( filename, 'r' );
data = textscan( fid, fmt, Inf, 'Delimiter', ',','headerlines',2) ;
fclose( fid ) ;

%% subset variables
% now using the _relArea.csv list (that includes relative area change)
varList = [8 12 13 28 38 51 53 54 39 46 59 45 4 29 7 11 64 26 2 3 21];

dArea = data{1}; % lake area change [m^2]
dAreaRelEnd = data{55};
dAreaRelStart = data{56};
detachInd = data{57};
attachInd = ~detachInd;

% in Hannah's Master_Lake_Data_Numerical_200717_withRelArea.csv
% 8 = temp change
% 12 = precip change
% 13 = initial lake area
% 28 = summed mass balance
% 38 = glacier response time
% 39 = mass blaance gradient
% 46 = mean near terminal thickness (I THINK)
% 51 = topgraphic setting (think 1-2 = A1-2, 3-4 = B1-2, 7 = E)
% 53 = initial lake area
% 54 = final area
% 59 = formed since record start? y/n

%varList = [3 4 5 8 9 12 13]; % variables to pull (column numbers)
%varList = [8 12 13 28 53 54];

% in Hannah's Master_Lake_Data_Numerical_200717.csv
% 8 = temp change
% 12 = precip change
% 13 = initial lake area
% 28 = summed mass balance
% 53 = initial lake area
% 54 = detached? 1 for yes
% 45 = median near terminal thickness
% 4 = elevation
% 29 = glacier area [km2]
% 7 = JJA temp 2000s
% 11 = DJF precip 2000s
% 64 = coast dist
% 26 = 2010s mass bal

predSub = []; % make empty storage matrix
for i = 1:length(varList) % iterate over all variables to pull
    
    numberNow = varList(i); % variable number (column number) to pull now
    varNow = data{numberNow}; % current varible
    
    if i == 1
        predSub = varNow; % create predSub (subset of predictor variables) if first iteration
    else
        predSub = [predSub varNow]; % otehrwise append to existing matrix
    end
end

%% rename variables

tempChange = predSub(:,1);
precipChange = predSub(:,2);
cumMassBal = predSub(:,4);
respTime = predSub(:,5);
%detachInd = predSub(:,6);
topoSet = predSub(:,6);
initArea = predSub(:,7);
finalArea = predSub(:,8);
dbdz = predSub(:,9);
hterm = predSub(:,10);
formedInd = predSub(:,11);
existedInd = ~formedInd;
medThick = predSub(:,12);
elev = predSub(:,13);
glArea = predSub(:,14);

temp = predSub(:,15);
precip = predSub(:,16);
coastDist = predSub(:,17);
b2010 = predSub(:,18);

long = predSub(:,19);
lat = predSub(:,20);
glWidth = predSub(:,21);



%% plot

figure(1)
clf
orient landscape

subplot(2,1,1)
hold on
posInd = dArea > 0;
negInd = dArea <= 0;
posAttachInd = posInd == 1 & ~detachInd == 1;
negAttachInd = negInd == 1 & ~detachInd == 1;


l = scatter(tempChange(posAttachInd),precipChange(posAttachInd),80,log10(dArea(posAttachInd)),'filled','linewidth',1);
l.MarkerEdgeColor = 'k';

m = scatter(tempChange(negAttachInd),precipChange(negAttachInd),100,log10(dArea(negAttachInd)),'linewidth',2);
m.Marker = 's';

c = colorbar();
ylabel(c,'log_{10}(|\Delta area|) [m^2]','fontsize',16)

yl = ylim();
xl = xlim();

plot(xl,[0 0],'k-')
set(gca,'fontsize',14)
grid on
xlabel('JJA temperature change [\circC]','fontsize',18)
ylabel('DJF precipitation change [mm]','fontsize',18)

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('tempPrecipChange_areaChange_withAbs.pdf','-dpdf','-r300')


%% mass balance vs. temp change

figure(2)
clf
orient landscape
hold on

scatter(tempChange(attachInd),precipChange(attachInd),80,cumMassBal(attachInd),'filled','linewidth',1);
colormap(flipud(parula))
c = colorbar();
ylabel(c,'Cum. mass balance [m w.e.]','fontsize',16)

yl = ylim();
xl = xlim();

plot(xl,[0 0],'k-')
set(gca,'fontsize',14)
grid on
xlabel('JJA temperature change [\circC]','fontsize',18)
ylabel('DJF precipitation change [mm]','fontsize',18)
title('climate params vs. cum mass bal')

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('tempPrecipChange_cumMassBal.pdf','-dpdf','-r300')



%% response time

figure(3)
clf
orient landscape
hold on

bInd = topoSet == 3 | topoSet ==4 & attachInd ==1;

l1 = scatter(respTime(attachInd),dArea(attachInd)/1e6,80,'b','linewidth',1);
hold on
l2 = scatter(respTime(bInd),dArea(bInd)/1e6,80,'r','filled','linewidth',1);
plot(xl,[0 0],'k-')

legend([l1,l2],{'All lakes','Type B'},'location','northwest')

yl = ylim();
xl = xlim();

set(gca,'fontsize',14)
grid on
xlabel('Glacier response time [a]','fontsize',18)
ylabel('Area change [km^2]','fontsize',18)

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('responseTimeAreaChange.pdf','-dpdf','-r300')


%% mass balance and relative area change

figure(4)
clf
orient landscape
subplot(2,1,1)
hold on

aInd = topoSet == 1 | topoSet ==2 & attachInd ==1;
eInd = topoSet == 7 & attachInd ==1;
proglInd = topoSet < 7 & attachInd ==1;


l1 = scatter(cumMassBal(proglInd),dAreaRelEnd(proglInd)*100,80,'k','linewidth',1);
l2 = scatter(cumMassBal(aInd),dAreaRelEnd(aInd)*100,80,'r','filled','linewidth',1);
l3 = scatter(cumMassBal(eInd),dAreaRelEnd(eInd)*100,80,'bs','filled','linewidth',1);
%l2 = scatter(respTime(bInd),dArea(bInd)/1e6,80,'r','filled','linewidth',1);
plot(xl,[0 0],'k-')
plot([0 0],yl,'k-')


legend([l1,l2,l3],{'Proglacial lakes','Proglacial, Type A','Ice-dammed, Type E'},'location','northeast')
ylim([-200 105])
yl = ylim();
xl = xlim();


set(gca,'fontsize',14)
grid on
xlabel('Cumulative mass balance [m w.e.]','fontsize',18)
ylabel('Relative area change, \Delta A / A_f [%]','fontsize',18)

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('relAreaCumMassBal.pdf','-dpdf','-r300')

%% mass balance gradient

figure(5)
clf
orient landscape
subplot(2,1,1)
hold on

aInd = topoSet == 1 | topoSet ==2 & attachInd ==1;
eInd = topoSet == 7 & attachInd ==1;
proglInd = topoSet < 7 & attachInd ==1;

l1 = scatter(dbdz(proglInd),dAreaRelEnd(proglInd),80,'k','linewidth',1);
l2 = scatter(dbdz(aInd),dAreaRelEnd(aInd),80,'r','filled','linewidth',1);
%l3 = scatter(cumMassBal(eInd),dAreaRelStart(eInd),80,'bs','filled','linewidth',1);
%l2 = scatter(respTime(bInd),dArea(bInd)/1e6,80,'r','filled','linewidth',1);

%plot([0 0],yl,'k-')


%legend([l1,l2,l3],{'Proglacial lakes','Progl., Type A','Ice-dammed'},'location','northeast')
ylim([-1.5 1.75])
xlim([0.25 1.5])
yl = ylim();
xl = xlim();

plot([0 2],[0 0],'k-')

set(gca,'fontsize',14)
grid on
xlabel('Balance gradient \times 100 [m w.e. m^{-1}]','fontsize',18)
ylabel('Relative area change, \Delta A / A_f [%]','fontsize',18)

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('relAreaMassBalGradient.pdf','-dpdf','-r300')


%% near terminal ice thickness

figure(5)
clf
orient landscape
subplot(2,1,1)
hold on

aInd = topoSet == 1 | topoSet ==2 & attachInd ==1;
eInd = topoSet == 7 & attachInd ==1;
proglInd = topoSet < 7 & attachInd ==1;

l1 = scatter(hterm(proglInd),dAreaRelEnd(proglInd)*100,80,'k','linewidth',1);
l2 = scatter(hterm(eInd),dAreaRelEnd(eInd)*100,80,'r','filled','linewidth',1);


%plot([0 0],yl,'k-')


legend([l1,l2],{'Proglacial lakes','Ice-dammed'},'location','southwest')
ylim([-250 105])
xlim([0 650])


set(gca,'fontsize',14)
grid on

ylabel('Relative area change, \Delta A / A_f [%]','fontsize',18)

subplot(2,1,2)
hold on

aInd = topoSet == 1 | topoSet ==2 & attachInd ==1;
eInd = topoSet == 7 & attachInd ==1;
proglInd = topoSet < 7 & attachInd ==1;

l1 = scatter(hterm(proglInd),dAreaRelEnd(proglInd)*100,80,'k','linewidth',1);
l2 = scatter(hterm(eInd),dAreaRelEnd(eInd)*100,80,'r','filled','linewidth',1);
xlim([0 650])

%plot([0 2],[0 0],'k-')

set(gca,'fontsize',14)
grid on
xlabel('Terminal ice thickness [m]','fontsize',18)
ylabel('Relative area change, \Delta A / A_f [%]','fontsize',18)

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('nearTermIceThick_relArea_exampleOutliersDrivingCorrelation.pdf','-dpdf','-r300')

%% mass balance and (end and init) relative area change

figure(4)
clf
orient landscape
subplot(2,1,1)
hold on

aInd = topoSet == 1 | topoSet ==2 & attachInd ==1;
eInd = topoSet == 7 & attachInd ==1;
proglInd = topoSet < 7 & attachInd ==1;
proglAndNewInd = proglInd == 1 & existedInd == 0;
proglAndExistInd = proglInd == 1 & existedInd == 1;
dammedAndNewInd = eInd == 1 & existedInd == 0;
dammedAndExistInd = eInd == 1 & existedInd == 1;

l1 = scatter(cumMassBal(proglInd),dAreaRelEnd(proglInd)*100,80,'k','linewidth',1);
l2 = scatter(cumMassBal(aInd),dAreaRelEnd(aInd)*100,80,'r','filled','linewidth',1);
l3 = scatter(cumMassBal(eInd),dAreaRelEnd(eInd)*100,80,'bs','filled','linewidth',1);
%l2 = scatter(respTime(bInd),dArea(bInd)/1e6,80,'r','filled','linewidth',1);
plot(xl,[0 0],'k-')
plot([0 0],yl,'k-')


legend([l1,l2,l3],{'Proglacial lakes','Proglacial, Type A','Ice-dammed, Type E'},'location','northeast')
ylim([-200 105])
yl = ylim();
xl = xlim();


set(gca,'fontsize',14)
grid on
ylabel('Relative area change, \Delta A / A_f [%]','fontsize',18)


subplot(2,1,2)
hold on
l1 = scatter(cumMassBal(proglInd),dAreaRelStart(proglInd)*100,80,'k','linewidth',1);
l2 = scatter(cumMassBal(aInd),dAreaRelStart(aInd)*100,80,'r','filled','linewidth',1);
l3 = scatter(cumMassBal(eInd),dAreaRelStart(eInd)*100,80,'bs','filled','linewidth',1);
%l2 = scatter(respTime(bInd),dArea(bInd)/1e6,80,'r','filled','linewidth',1);
plot(xl,[0 0],'k-')
plot([0 0],yl,'k-')


legend([l1,l2,l3],{'Proglacial lakes','Proglacial, Type A','Ice-dammed, Type E'},'location','northeast')
ylim([-105 200])
yl = ylim();
xl = xlim();


set(gca,'fontsize',14)
grid on
xlabel('Cumulative mass balance [m w.e.]','fontsize',18)
ylabel('Relative area change, \Delta A / A_i [%]','fontsize',18)

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('relAreaCumMassBal.pdf','-dpdf','-r300')

%% start vs end area plot with relative change lines

figure(8)
clf
orient landscape
hold on

l1 = plot(initArea(proglAndExistInd)/1e6,finalArea(proglAndExistInd)/1e6,'b.','markersize',25);
scatter(initArea(proglAndNewInd)/1e6,finalArea(proglAndNewInd)/1e6,40,'bo');
l1 = scatter(initArea(proglAndExistInd)/1e6,finalArea(proglAndExistInd)/1e6,40,'bo','filled');
l2 = scatter(initArea(dammedAndExistInd)/1e6,finalArea(dammedAndExistInd)/1e6,40,'rd','filled');
scatter(initArea(dammedAndNewInd)/1e6,finalArea(dammedAndNewInd)/1e6,40,'rd');
plot([0 100],[0 100],'k--','linewidth',2)
plot([0 100],[0 50],'k:','linewidth',1)
plot([0 100],[0 200],'k:','linewidth',1)
plot([0 100],[0 150],'k:','linewidth',1)
plot([0 100],[0 75],'k:','linewidth',1)

areaThresh = 0.0001 * 1e6; % area threshold (m2) for defining a lake as existing or not (zero areas actually have some area b/c need to click a few times to define observation as zero)
existFromStartInd = initArea > areaThresh;
sum(existFromStartInd)

%legend([l1,l2],{'Proglacial','Ice-dammed'},'location','northwest')
xlabel('Starting area [km^2]','fontsize',16)
ylabel('Ending area [km^2]','fontsize',16)
grid on
set(gca,'fontsize',14)
% zoom out
if 1
ylim([0 100])
xlim([0 60])
end

% zoom in
if 0
xlim([0 10])
ylim([0 20])
end

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('startAreaEndArea_withGuideLines_zoomOut_withNewLakes.pdf','-dpdf','-r300')

%% start vs end area plot with relative change lines - log scale

xl = [1e-2 1e2];
yl = [1e-2 1e2];


figure(88)
clf
orient landscape
hold on

l1 = plot(initArea(proglAndExistInd)/1e6,finalArea(proglAndExistInd)/1e6,'b.','markersize',25);
scatter(initArea(proglAndNewInd)/1e6,finalArea(proglAndNewInd)/1e6,40,'bo');
l1 = scatter(initArea(proglAndExistInd)/1e6,finalArea(proglAndExistInd)/1e6,40,'bo','filled');
l2 = scatter(initArea(dammedAndExistInd)/1e6,finalArea(dammedAndExistInd)/1e6,40,'rd','filled');
scatter(initArea(dammedAndNewInd)/1e6,finalArea(dammedAndNewInd)/1e6,40,'rd');
plot([1e-6 100],[1e-6 100],'k--','linewidth',2)
plot([1e-6 100],[1e-6 50],'k:','linewidth',1)
plot([1e-6 100],[1e-6 200],'k:','linewidth',1)
plot([1e-6 100],[1e-6 150],'k:','linewidth',1)
plot([1e-6 100],[1e-6 75],'k:','linewidth',1)
set(gca,'xscale','log')
set(gca,'yscale','log')
%xlim(xl)
%ylim(yl)



%legend([l1,l2],{'Proglacial','Ice-dammed'},'location','northwest')
xlabel('Starting area [km^2]','fontsize',16)
ylabel('Ending area [km^2]','fontsize',16)
grid on
set(gca,'fontsize',14)


h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print('startAreaEndArea_withGuideLines_zoomOut_withNewLakes_logLog.pdf','-dpdf','-r300')



%% test for how correlation depends on relative area metric used, statistical test used
predictorNow = cumMassBal;

predInd =  ~isnan(predictorNow);
predAndExistInd = predInd == 1 & existedInd ==1;
predAndNewInd = predInd == 1 & existedInd == 0;
predAndProglInd = predInd == 1 & proglInd == 1;
predAndDammedInd = predInd == 1 & eInd == 1;

predLims = [nanmin(predictorNow),nanmax(predictorNow)];

pFitRelEndWithAll = polyfit(cumMassBal(predInd),dAreaRelEnd(predInd),1);
pFitRelEndExisted = polyfit(cumMassBal(predAndExistInd),dAreaRelEnd(predAndExistInd),1);
pFitRelStartWithAll = polyfit(cumMassBal(predInd),dAreaRelStart(predInd),1);
pFitRelStartExisted = polyfit(cumMassBal(predAndExistInd),dAreaRelStart(predAndExistInd),1);
pFitRelEndProgl = polyfit(cumMassBal(predAndProglInd),dAreaRelEnd(predAndProglInd),1);
pFitRelEndDammed = polyfit(cumMassBal(predAndDammedInd),dAreaRelEnd(predAndDammedInd),1);

fitLineRelEndWithAll = polyval(pFitRelEndWithAll,predLims);
fitLineRelEndExisted = polyval(pFitRelEndExisted,predLims);
fitLineRelStartWithAll = polyval(pFitRelStartWithAll,predLims);
fitLineRelStartExisted = polyval(pFitRelStartExisted,predLims);
fitLineRelEndProgl = polyval(pFitRelEndProgl,predLims);
fitLineRelEndDammed = polyval(pFitRelEndDammed,predLims);

[mFitRelStartAll, bFitRelStartAll] = TheilSen([predictorNow(predInd),dAreaRelStart(predInd)]);
senFitRelStartAll = predLims*mFitRelStartAll + bFitRelStartAll;

[mFitRelEndAll, bFitRelEndAll] = TheilSen([predictorNow(predInd),dAreaRelEnd(predInd)]);
senFitRelEndAll = predLims*mFitRelEndAll + bFitRelEndAll;

[mFitRelEndProgl, bFitRelEndProgl] = TheilSen([predictorNow(predAndProglInd),dAreaRelEnd(predAndProglInd)]);
senFitRelEndProgl = predLims*mFitRelEndProgl + bFitRelEndProgl;

[mFitRelEndDammed, bFitRelEndDammed] = TheilSen([predictorNow(predAndDammedInd),dAreaRelEnd(predAndDammedInd)]);
senFitRelEndDammed = predLims*mFitRelEndDammed + bFitRelEndDammed;


%[rEndAll,rpEndAll] = corr(cumMassBal(predInd),dAreaRelEnd(predInd),'type','pearson');
%[tauEndAll,taupEndAll] = corr(cumMassBal(predInd),dAreaRelEnd(predInd),'type','kendall');
%[rStartAll,rpStartAll] = corr(cumMassBal(predInd),dAreaRelStart(predInd),'type','pearson');
%[tauStartAll,taupStartAll] = corr(cumMassBal(predInd),dAreaRelStart(predInd),'type','kendall');

%[rEndProgl,rpEndProgl] = corr(cumMassBal(predAndProglInd),dAreaRelEnd(predAndProglInd),'type','pearson');
%[tauEndProgl,taupEndProgl] = corr(cumMassBal(predAndProglInd),dAreaRelEnd(predAndProglInd),'type','kendall');

%[rEndDammed,rpEndDammed] = corr(cumMassBal(predAndDammedInd),dAreaRelEnd(predAndDammedInd),'type','pearson');
%[tauEndDammed,taupEndDammed] = corr(cumMassBal(predAndDammedInd),dAreaRelEnd(predAndDammedInd),'type','kendall');


figure(13)
clf
orient landscape
hold on

h1 = plot(cumMassBal(predAndProglInd),dAreaRelEnd(predAndProglInd),'b.','markersize',20);
h2 = plot(cumMassBal(predAndDammedInd),dAreaRelEnd(predAndDammedInd),'rd','markerfacecolor','r','markersize',8);
h3 = plot(predLims,fitLineRelEndProgl,'b--');
h4 = plot(predLims,senFitRelEndProgl,'b');
h5 = plot(predLims,fitLineRelEndDammed,'r--');
h6 = plot(predLims,senFitRelEndDammed,'r');
ylim([-5 2])
xlim(predLims)
legend([h1,h2,h3,h4],{'Progl.','Ice-dammed','polyfit','theil-sen'},'location','southeast')
xlabel('Cum. mass bal. [m w.e.]','fontsize',16)
ylabel('Relative area change, \DeltaA/A_f [-]','fontsize',16)
grid on


h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('cumulativeMassBal_relAreaChagne_iceDammed_proglacial_polyfitAndSenSlope.pdf','-dpdf','-r300')

%%

figure(12)
clf
orient landscape

subplot(2,1,1)
hold on
h1 = plot(cumMassBal(predInd),dAreaRelEnd(predInd),'r.','markersize',20);
h2 = plot(cumMassBal(predAndExistInd),dAreaRelEnd(predAndExistInd),'b.','markersize',20);
plot(predLims,fitLineRelEndWithAll,'r')
plot(predLims,fitLineRelEndExisted,'b')
plot(predLims,senFitRelEndAll,'k')
legend([h1,h2],{'All lakes','Existing lakes'},'location','southeast')
ylabel('\DeltaA/A_f','fontsize',16)
set(gca,'fontsize',14)
grid on
ylim([-5 1.5])
xlim(predLims)


subplot(2,1,2)
hold on
h1 = plot(cumMassBal(predInd),dAreaRelStart(predInd),'r.','markersize',20);
h2 = plot(cumMassBal(predAndExistInd),dAreaRelStart(predAndExistInd),'b.','markersize',20);
plot(predLims,fitLineRelStartWithAll,'r')
plot(predLims,fitLineRelStartExisted,'b')
plot(predLims,senFitRelStartAll,'k')
ylim([-1.5 10])
xlim(predLims)

set(gca,'fontsize',14)
ylabel('\DeltaA/A_i','fontsize',16)
xlabel('Cum. mass bal. [m w.e.]','fontsize',16)
grid on

%% recreating some of Hannah's correlations to see effect of non paramteric stats


predictorNow = cumMassBal;

predInd =  ~isnan(predictorNow);
predAndExistInd = predInd == 1 & existedInd ==1;
predAndProglExistInd = predInd == 1 & existedInd ==1 & proglInd == 1;
predAndProglNewInd = predInd == 1 & existedInd ==0 & proglInd == 1;
predAndProglInd = predInd == 1 & proglInd == 1;
%predAndDamemdExistInd = predInd == 1 & existedInd ==1 & eInd == 1;
predAndDammedExistInd = predInd == 1 & existedInd ==1 & eInd == 1;

predAndDammedNewInd = predInd == 1 & existedInd == 0 & eInd == 1;

%predLims = [nanmin(predictorNow),nanmax(predictorNow)];


%% relative area vs starting area

figure(11)
clf
orient landscape
% hold on

absLim = [- 15 30];
relLim = [-100 500];
xl = [0 55];
marksize = 60;

[mInitAbsP, bInitAbsP] = TheilSen([initArea(predAndProglExistInd)/1e6,dArea(predAndProglExistInd)/1e6]);
[mInitAbsD, bInitAbsD] = TheilSen([initArea(predAndDammedExistInd)/1e6,dArea(predAndDammedExistInd)/1e6]);
[mInitRelP, bInitRelP] = TheilSen([initArea(predAndProglExistInd)/1e6,dAreaRelStart(predAndProglExistInd)*100]);
[mInitRelD, bInitRelD] = TheilSen([initArea(predAndDammedExistInd)/1e6,dAreaRelStart(predAndDammedExistInd)*100]);


subplot(2,1,1)
hold on
scatter(initArea(predAndProglExistInd)/1e6, dArea(predAndProglExistInd)/1e6,marksize,'bo','filled') %corr1
scatter(initArea(predAndProglNewInd)/1e6, dArea(predAndProglNewInd)/1e6,marksize,'bo') %corr1
scatter(initArea(predAndDammedExistInd)/1e6, dArea(predAndDammedExistInd)/1e6,marksize,'rd','filled')
scatter(initArea(predAndDammedNewInd)/1e6, dArea(predAndDammedNewInd)/1e6,marksize,'rd')
plot(xl,[0 0],'k:','linewidth',2)
plot(xl,xl.*mInitAbsP + bInitAbsP,'b','linewidth',2)
plot(xl,xl.*mInitAbsD + bInitAbsD,'r','linewidth',2)

grid on
xlim(xl)
%ylim(absLim)
set(gca,'fontsize',14)
ylabel('Absolute area change, \DeltaA [km^2]','fontsize',16)
%legend('Proglacial','Ice-dammed','location','southeast')

subplot(2,1,2)
hold on
scatter(initArea(predAndProglExistInd)/1e6, dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
%scatter(initArea(predAndProglNewInd)/1e6, dAreaRelStart(predAndProglNewInd)*100,marksize,'bo') %corr1
scatter(initArea(predAndDammedExistInd)/1e6, dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
z = plot(xl,[0 0],'k:','linewidth',2);
t1 = plot(xl,xl.*mInitRelP + bInitRelP,'b','linewidth',2);
t2 = plot(xl,xl.*mInitRelD + bInitRelD,'r','linewidth',2);

%legend([t1,t2,z],{'Theil-Sen: proglacial','Theil-Sen: ice-dammed','No change'},'location','northeast');

grid on
xlim(xl)
%ylim(relLim)
%semilogy(initArea(corrInd)/1e6, dAreaRelStart(corrInd)*-100,'b.','markersize',20) %corr1
%plot(predictorNow(existAndCorrInd), dAreaRelStart(existAndCorrInd),'b.','markersize',20) %corr1
%title('A_i & all')
xlabel('Starting area [km^2]','fontsize',16)
ylabel('Relative area change, \DeltaA/A_i [%]','fontsize',16)
set(gca,'fontsize',14)
grid on

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('intArea_AbsAndRelChange_showNewLakes.pdf','-dpdf','-r300')


%% median near terminal ice thickenss


[mMedThickAbsP, bMedThickAbsP] = TheilSen([medThick(predAndProglExistInd),dArea(predAndProglExistInd)/1e6]);
[mMedThickAbsD, bMedThickAbsD] = TheilSen([medThick(predAndDammedExistInd),dArea(predAndDammedExistInd)/1e6]);
[mMedThickRelP, bMedThickRelP] = TheilSen([medThick(predAndProglExistInd),dAreaRelStart(predAndProglExistInd)*100]);
[mMedThickRelD, bMedThickRelD] = TheilSen([medThick(predAndDammedExistInd),dAreaRelStart(predAndDammedExistInd)*100]);

xl = [40 700];
%yl1 = [-15 50];
%yl2 = [-100 800];

figure(30)
clf
orient landscape

subplot(2,1,1)
hold on
scatter(medThick(predAndProglExistInd), dArea(predAndProglExistInd)/1e6,marksize,'bo','filled') %corr1
scatter(medThick(predAndProglNewInd), dArea(predAndProglNewInd)/1e6,marksize,'bo') %corr1
scatter(medThick(predAndDammedExistInd), dArea(predAndDammedExistInd)/1e6,marksize,'rd','filled')
scatter(medThick(predAndDammedNewInd), dArea(predAndDammedNewInd)/1e6,marksize,'rd','filled')
plot(xl,[0 0],'k:','linewidth',2)
t1 = plot(xl,xl.*mMedThickAbsP + bMedThickAbsP,'b','linewidth',2);
plot(xl,xl.*mMedThickAbsD + bMedThickAbsD,'r--','linewidth',0.5)

grid on
xlim(xl)
%ylim(absLim)
set(gca,'fontsize',14)
ylabel('Absolute area change, \DeltaA [km^2]','fontsize',16)
%legend('Proglacial','Ice-dammed','location','northeast')

subplot(2,1,2)
hold on
scatter(medThick(predAndProglExistInd), dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
scatter(medThick(predAndProglNewInd), dAreaRelStart(predAndProglNewInd)*100,marksize,'bo') %corr1
scatter(medThick(predAndDammedExistInd), dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
plot(xl,[0 0],'k:','linel',2)
t2 = plot(xl,xl.*mMedThickRelP + bMedThickRelP,'b--','linewidth',0.5);
t3 = plot(xl,xl.*mMedThickRelD + bMedThickRelD,'r','linewidth',0.5);

%legend([t2,t3,t1],{'p > 0.1','0.05 < p \leq 0.10','p\leq 0.05'},'location','northeast')
grid on
xlim(xl)
%ylim(relLim)
set(gca,'fontsize',14)
xlabel('Median lake-adjacent glacier thickness [m]','fontsize',16)
ylabel('Relative area change, \DeltaA/A_i [%]','fontsize',16)
set(gca,'fontsize',14)
grid on

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('medIceThick_AbsAndRelChange_withNewLakes.pdf','-dpdf','-r300')

%% elevation


[mElevAbsP, bElevAbsP] = TheilSen([elev(predAndProglExistInd),dArea(predAndProglExistInd)/1e6]);
[mElevAbsD, bElevAbsD] = TheilSen([elev(predAndDammedExistInd),dArea(predAndDammedExistInd)/1e6]);
[mElevRelP, bElevRelP] = TheilSen([elev(predAndProglExistInd),dAreaRelStart(predAndProglExistInd)*100]);
[mElevRelD, bElevRelD] = TheilSen([elev(predAndDammedExistInd),dAreaRelStart(predAndDammedExistInd)*100]);

xl = [0 2200];


figure(31)
clf
orient landscape

subplot(2,1,1)
hold on
scatter(elev(predAndProglExistInd), dArea(predAndProglExistInd)/1e6,marksize,'bo','filled') %corr1
scatter(elev(predAndProglNewInd), dArea(predAndProglNewInd)/1e6,marksize,'bo') %corr1
scatter(elev(predAndDammedExistInd), dArea(predAndDammedExistInd)/1e6,marksize,'rd','filled')
plot(xl,[0 0],'k:','linewidth',2)
plot(xl,xl.*mElevAbsP + bElevAbsP,'b','linewidth',2);
plot(xl,xl.*mElevAbsD + bElevAbsD,'r--','linewidth',0.5)

grid on
xlim(xl)
ylim(absLim)
set(gca,'fontsize',14)
ylabel('Absolute area change, \DeltaA [km^2]','fontsize',16)
%legend('Proglacial','Ice-dammed','location','northeast')

subplot(2,1,2)
hold on
scatter(elev(predAndProglExistInd), dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
scatter(elev(predAndProglNewInd), dAreaRelStart(predAndProglNewInd)*100,marksize,'bo')
scatter(elev(predAndDammedExistInd), dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
plot(xl,[0 0],'k:','linewidth',2)
plot(xl,xl.*mElevRelP + bElevRelP,'b-','linewidth',2);
plot(xl,xl.*mElevRelD + bElevRelD,'r--','linewidth',0.5);

%legend([t2,t3,t1],{'p > 0.1','0.05 < p \leq 0.10','p\leq 0.05'},'location','northeast')
grid on
xlim(xl)
ylim(relLim)
set(gca,'fontsize',14)
xlabel('Lake elevation [m a.s.l.]','fontsize',16)
ylabel('Relative area change, \DeltaA/A_i [%]','fontsize',16)
set(gca,'fontsize',14)
grid on

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('elev_AbsAndRelChange_wideLims.pdf','-dpdf','-r300')


%% glacier area
[mAreaAbsP, bAreaAbsP] = TheilSen([glArea(predAndProglExistInd),dArea(predAndProglExistInd)/1e6]);
[mAreaAbsD, bAreaAbsD] = TheilSen([glArea(predAndDammedExistInd),dArea(predAndDammedExistInd)/1e6]);
[mAreaRelP, bAreaRelP] = TheilSen([glArea(predAndProglExistInd),dAreaRelStart(predAndProglExistInd)*100]);
[mAreaRelD, bAreaRelD] = TheilSen([glArea(predAndDammedExistInd),dAreaRelStart(predAndDammedExistInd)*100]);

xl = [0 1100];


figure(31)
clf
orient landscape

subplot(2,1,1)
hold on
scatter(glArea(predAndProglExistInd), dArea(predAndProglExistInd)/1e6,marksize,'bo','filled') %corr1
scatter(glArea(predAndProglNewInd), dArea(predAndProglNewInd)/1e6,marksize,'bo') %corr1
scatter(glArea(predAndDammedExistInd), dArea(predAndDammedExistInd)/1e6,marksize,'rd','filled')
plot(xl,[0 0],'k:','linewidth',2)
plot(xl,xl.*mAreaAbsP + bAreaAbsP,'b','linewidth',2);
plot(xl,xl.*mAreaAbsD + bAreaAbsD,'r--','linewidth',0.5)

grid on
xlim(xl)
ylim(absLim)
set(gca,'fontsize',14)
ylabel('Absolute area change, \DeltaA [km^2]','fontsize',16)
%legend('Proglacial','Ice-dammed','location','northeast')

subplot(2,1,2)
hold on
scatter(glArea(predAndProglExistInd), dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
scatter(glArea(predAndProglNewInd), dAreaRelStart(predAndProglNewInd)*100,marksize,'bo')
scatter(glArea(predAndDammedExistInd), dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
plot(xl,[0 0],'k:','linewidth',2)
plot(xl,xl.*mAreaRelP + bAreaRelP,'b--','linewidth',0.5);
plot(xl,xl.*mAreaRelD + bAreaRelD,'r--','linewidth',0.5);

%legend([t2,t3,t1],{'p > 0.1','0.05 < p \leq 0.10','p\leq 0.05'},'location','northeast')
grid on
xlim(xl)
ylim(relLim)
set(gca,'fontsize',14)
xlabel('Glacier area [km^2]','fontsize',16)
ylabel('Relative area change, \DeltaA/A_i [%]','fontsize',16)
set(gca,'fontsize',14)
grid on

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('glArea_AbsAndRelChange.pdf','-dpdf','-r300')



%% jja temp  and djf precip

xl1 = [5 15];
xl2 = [0 2];
xl3 = [0 1200];
xl4 = [-300 250];

figure(33)
clf
orient landscape

subplot(4,1,1)
hold on
scatter(temp(predAndProglExistInd), dArea(predAndProglExistInd)/1e6,marksize,'bo','filled') %corr1
scatter(temp(predAndProglNewInd), dArea(predAndProglNewInd)/1e6,marksize,'bo') %corr1
scatter(temp(predAndDammedExistInd), dArea(predAndDammedExistInd)/1e6,marksize,'rd','filled')
plot(xl1,[0 0],'k:','linewidth',2)
%plot(xl,xl.*mAreaAbsP + bAreaAbsP,'b','linewidth',2);
%plot(xl,xl.*mAreaAbsD + bAreaAbsD,'r--','linewidth',0.5)
xlabel('JJA temperature [^{\circ}C]')
ylabel('\DeltaA [km^2]')
ylim(absLim)
grid on

subplot(4,1,2)
hold on
scatter(tempChange(predAndProglExistInd), dArea(predAndProglExistInd)/1e6,marksize,'bo','filled') %corr1
scatter(tempChange(predAndProglNewInd), dArea(predAndProglNewInd)/1e6,marksize,'bo') %corr1
scatter(tempChange(predAndDammedExistInd), dArea(predAndDammedExistInd)/1e6,marksize,'rd','filled')
plot(xl2,[0 0],'k:','linewidth',2)
%plot(xl,xl.*mAreaAbsP + bAreaAbsP,'b','linewidth',2);
%plot(xl,xl.*mAreaAbsD + bAreaAbsD,'r--','linewidth',0.5)
xlabel('JJA temperature change [^{\circ}C]')
ylim(absLim)
ylabel('\DeltaA [km^2]')
grid on


subplot(4,1,3)
hold on
scatter(precip(predAndProglExistInd), dArea(predAndProglExistInd)/1e6,marksize,'bo','filled') %corr1
scatter(precip(predAndProglNewInd), dArea(predAndProglNewInd)/1e6,marksize,'bo') %corr1
scatter(precip(predAndDammedExistInd), dArea(predAndDammedExistInd)/1e6,marksize,'rd','filled')
plot(xl3,[0 0],'k:','linewidth',2)
%plot(xl,xl.*mAreaAbsP + bAreaAbsP,'b','linewidth',2);
%plot(xl,xl.*mAreaAbsD + bAreaAbsD,'r--','linewidth',0.5)
xlabel('DJF precipitation [mm w.e.]')
ylabel('\DeltaA [km^2]')
ylim(absLim)
xlim(xl3)
grid on

subplot(4,1,4)
hold on
scatter(precipChange(predAndProglExistInd), dArea(predAndProglExistInd)/1e6,marksize,'bo','filled') %corr1
scatter(precipChange(predAndProglNewInd), dArea(predAndProglNewInd)/1e6,marksize,'bo') %corr1
scatter(precipChange(predAndDammedExistInd), dArea(predAndDammedExistInd)/1e6,marksize,'rd','filled')
plot(xl4,[0 0],'k:','linewidth',2)
%plot(xl,xl.*mAreaAbsP + bAreaAbsP,'b','linewidth',2);
%plot(xl,xl.*mAreaAbsD + bAreaAbsD,'r--','linewidth',0.5)

xlabel('DJF precipitation change [mm w.e.]')
ylim(absLim)
xlim(xl4)
ylabel('\DeltaA [km^2]')
grid on


h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print('tempPrecipAndChagne_AbsChange.pdf','-dpdf','-r300')


%% climate params with relative area change
xl1 = [5 15];
xl2 = [0 2];
xl3 = [0 1200];
xl4 = [-300 250];

figure(35)
clf
orient landscape

subplot(4,1,1)
hold on
scatter(temp(predAndProglExistInd), dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
scatter(temp(predAndProglNewInd), dAreaRelStart(predAndProglNewInd)*100,marksize,'bo') %corr1
scatter(temp(predAndDammedExistInd), dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
plot(xl1,[0 0],'k:','linewidth',2)
%plot(xl,xl.*mAreaAbsP + bAreaAbsP,'b','linewidth',2);
%plot(xl,xl.*mAreaAbsD + bAreaAbsD,'r--','linewidth',0.5)
xlabel('JJA temperature [^{\circ}C]')
ylabel('\DeltaA/A_i [%]')
ylim(relLim)
grid on

subplot(4,1,2)
hold on
scatter(tempChange(predAndProglExistInd), dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
scatter(tempChange(predAndProglNewInd), dAreaRelStart(predAndProglNewInd)*100,marksize,'bo') %corr1
scatter(tempChange(predAndDammedExistInd), dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
plot(xl2,[0 0],'k:','linewidth',2)
%plot(xl,xl.*mAreaAbsP + bAreaAbsP,'b','linewidth',2);
%plot(xl,xl.*mAreaAbsD + bAreaAbsD,'r--','linewidth',0.5)
xlabel('JJA temperature change [^{\circ}C]')
ylim(relLim)
ylabel('\DeltaA/A_i [%]')
grid on


subplot(4,1,3)
hold on
scatter(precip(predAndProglExistInd), dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
scatter(precip(predAndProglNewInd), dAreaRelStart(predAndProglNewInd)*100,marksize,'bo') %corr1
scatter(precip(predAndDammedExistInd), dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
plot(xl3,[0 0],'k:','linewidth',2)
%plot(xl,xl.*mAreaAbsP + bAreaAbsP,'b','linewidth',2);
%plot(xl,xl.*mAreaAbsD + bAreaAbsD,'r--','linewidth',0.5)
xlabel('DJF precipitation [mm w.e.]')
ylabel('\DeltaA/A_i [%]')
ylim(relLim)
xlim(xl3)
grid on

subplot(4,1,4)
hold on
scatter(precipChange(predAndProglExistInd), dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
scatter(precipChange(predAndProglNewInd), dAreaRelStart(predAndProglNewInd)*100,marksize,'bo') %corr1
scatter(precipChange(predAndDammedExistInd), dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
plot(xl4,[0 0],'k:','linewidth',2)
%plot(xl,xl.*mAreaAbsP + bAreaAbsP,'b','linewidth',2);
%plot(xl,xl.*mAreaAbsD + bAreaAbsD,'r--','linewidth',0.5)

xlabel('DJF precipitation change [mm w.e.]')
ylim(relLim)
xlim(xl4)
ylabel('\DeltaA/A_i [%]')
grid on


h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('tempPrecipAndChagne_RelChange.pdf','-dpdf','-r300')

%% cum mass bal
xl1 = [5 15];
xl2 = [0 2];
xl3 = [0 1200];
xl4 = [-300 250];

figure(35)
clf
orient landscape

subplot(2,1,2)
hold on
scatter(cumMassBal(predAndProglExistInd), dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
scatter(cumMassBal(predAndProglNewInd), dAreaRelStart(predAndProglNewInd)*100,marksize,'bo') %corr1
scatter(cumMassBal(predAndDammedExistInd), dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
plot(xl1,[0 0],'k:','linewidth',2)
%plot(xl,xl.*mAreaAbsP + bAreaAbsP,'b','linewidth',2);
%plot(xl,xl.*mAreaAbsD + bAreaAbsD,'r--','linewidth',0.5)
xlabel('JJA temperature [^{\circ}C]')
ylabel('\DeltaA/A_i [%]')
ylim(relLim)
grid on

%% coast dist
xl = [0 700];

[tAbsPro, pAbsPro] = corr(coastDist(predAndProglExistInd),dArea(predAndProglExistInd)/1e6,'type','kendall');
[tRelPro, pRelPro] = corr(coastDist(predAndProglExistInd),dAreaRelStart(predAndProglExistInd)*1e2,'type','kendall');
[tAbsDam, pAbsDam] = corr(coastDist(predAndDammedExistInd),dArea(predAndDammedExistInd)/1e6,'type','kendall');
[tRelDam, pRelDam] = corr(coastDist(predAndDammedExistInd),dAreaRelStart(predAndDammedExistInd)*1e2,'type','kendall');

[mDistAbsP, bDistAbsP] = TheilSen([coastDist(predAndProglExistInd),dArea(predAndProglExistInd)/1e6]);
[mDistAbsD, bDistAbsD] = TheilSen([coastDist(predAndDammedExistInd),dArea(predAndDammedExistInd)/1e6]);
[mDistRelP, bDistRelP] = TheilSen([coastDist(predAndProglExistInd),dAreaRelStart(predAndProglExistInd)*100]);
[mDistRelD, bDistRelD] = TheilSen([coastDist(predAndDammedExistInd),dAreaRelStart(predAndDammedExistInd)*100]);


figure(36)
clf
orient landscape
subplot(2,1,1)
hold on
scatter(coastDist(predAndProglExistInd), dArea(predAndProglExistInd)/1e6,marksize,'bo','filled') %corr1
scatter(coastDist(predAndProglNewInd), dArea(predAndProglNewInd)/1e6,marksize,'bo') %corr1
scatter(coastDist(predAndDammedExistInd), dArea(predAndDammedExistInd)/1e6,marksize,'rd','filled')
plot(xl,[0 0],'k:','linewidth',2)
plot(xl,xl.*mDistAbsP + bDistAbsP,'b','linewidth',2);
plot(xl,xl.*mDistAbsD + bDistAbsD,'r--','linewidth',0.5)
%ylim(absLim)
grid on
set(gca,'fontsize',14)
ylabel('Absolute area change, \DeltaA [km^2]','fontsize',16)
xlim(xl)


subplot(2,1,2)
hold on
scatter(coastDist(predAndProglExistInd), dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
scatter(coastDist(predAndProglNewInd), dAreaRelStart(predAndProglNewInd)*100,marksize,'bo') %corr1
scatter(coastDist(predAndDammedExistInd), dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
plot(xl,[0 0],'k:','linewidth',2)
plot(xl,xl.*mDistRelP + bDistRelP,'b','linewidth',.05);
plot(xl,xl.*mDistRelD + bDistRelD,'r--','linewidth',0.5)

xlabel('Distance to coast [km]')
ylabel('Relative area change, \DeltaA/A_i [%]','fontsize',16)
%ylim(relLim)
grid on
set(gca,'fontsize',14)
xlim(xl)


h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('distToCoast_AbsAndRelChange.pdf','-dpdf','-r300')
%print('distToCoast_AbsAndRelChange.png','-dpng','-r300')

%% dist to coast and climate

figure(37)
clf
orient landscape

ax(1) = subplot(3,1,1);
scatter(coastDist,temp,30,tempChange,'filled','linewidth',0.5,'markeredgecolor','k')
caxis([0 1.5])
colormap(ax(1),hot)
cb1 = colorbar();
ylabel(cb1,'\DeltaTemp. [^{\circ}C]','fontsize',14)
ylabel('JJA temperature [^{\circ}C]','fontsize',16)
grid on
set(gca,'fontsize',14)

subplot(3,1,2)
scatter(coastDist,precip,30,precipChange,'filled','linewidth',0.5,'markeredgecolor','k')
caxis([-200 200])
cm2 = flipud(customcolormap_preset('brown-white-pool'));
colormap(cm2)
cb2 = colorbar();
ylabel(cb2,'\DeltaPrecip. [mm]','fontsize',14)

ylabel('DJF precipitation [mm]','fontsize',16)
grid on
set(gca,'fontsize',14)

ax(3) = subplot(3,1,3);
scatter(coastDist,elev,30,log10(abs(dArea)),'filled','linewidth',0.5,'markeredgecolor','k')
caxis([4 7])
colormap(ax(3),parula)
cb3 = colorbar();
ylabel(cb3,'log_{10}(|\DeltaA|) [m]','fontsize',14)
ylabel('Elevation [m a.s.l.]','fontsize',16)
grid on
set(gca,'fontsize',14)
xlabel('Distance to coast [km]')

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print('distToCoast_tempPrecipAreaAndChange.pdf','-dpdf','-r300')


%%
%% mass bal & grad
xl = [-2 1];
xl2 = [0.3 1.5];

[mBAbsP, bBAbsP] = TheilSen([b2010(predAndProglExistInd),dArea(predAndProglExistInd)/1e6]);
[mBAbsD, bBAbsD] = TheilSen([b2010(predAndDammedExistInd),dArea(predAndDammedExistInd)/1e6]);

[mDbdzRelP, bDbdzRelP] = TheilSen([dbdz(predAndProglExistInd),dAreaRelStart(predAndProglExistInd)*100]);
[mDbdzRelD, bDbdzRelD] = TheilSen([dbdz(predAndDammedExistInd),dAreaRelStart(predAndDammedExistInd)*100]);


figure(39)
clf
orient landscape
subplot(2,1,1)
hold on
scatter(b2010(predAndProglExistInd), dArea(predAndProglExistInd)/1e6,marksize,'bo','filled') %corr1
scatter(b2010(predAndProglNewInd), dArea(predAndProglNewInd)/1e6,marksize,'bo') %corr1
scatter(b2010(predAndDammedExistInd), dArea(predAndDammedExistInd)/1e6,marksize,'rd','filled')
plot(xl,[0 0],'k:','linewidth',2)
plot([0 0],absLim,'k-','linewidth',1)
plot(xl,xl.*mBAbsP + bBAbsP,'b','linewidth',2);
plot(xl,xl.*mBAbsD + bBAbsD,'r--','linewidth',0.5)
%ylim(absLim)
grid on
set(gca,'fontsize',14)
ylabel('Absolute area change, \DeltaA [km^2]','fontsize',16)
xlabel('2010s glacier mass balance [m w.e. a^{-1}]','fontsize',16)
xlim(xl)


subplot(2,1,2)
hold on
scatter(dbdz(predAndProglExistInd), dAreaRelStart(predAndProglExistInd)*100,marksize,'bo','filled') %corr1
scatter(dbdz(predAndProglNewInd), dAreaRelStart(predAndProglNewInd)*100,marksize,'bo') %corr1
scatter(dbdz(predAndDammedExistInd), dAreaRelStart(predAndDammedExistInd)*100,marksize,'rd','filled')
plot(xl2,[0 0],'k:','linewidth',2)
plot(xl2,xl2.*mDbdzRelP + bDbdzRelP,'b','linewidth',2);
plot(xl2,xl2.*mDbdzRelD + bDbdzRelD,'r--','linewidth',0.5)

xlabel('Mass balance gradient [m w.e. a^{-1} 100 m^{-1}]')
ylabel('Relative area change, \DeltaA/A_i [%]','fontsize',16)
%ylim(relLim)
grid on
set(gca,'fontsize',14)
xlim(xl2)


h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('massBalAndGradient_AbsAndRelChange.pdf','-dpdf','-r300')
%print('massBalAndGradient_AbsAndRelChange.png','-dpng','-r300')






%% principal components
% the /1e6 kind of things no longer needed - the lines below normalize the
% data. Before, just got hooked onto whatever variable had the biggest
% values

predSubCleaned = [lat, long, elev, initArea/1e6, coastDist, temp, tempChange, precip, precipChange, glArea, glWidth/1e3, hterm, dbdz, b2010, cumMassBal ];
%                               1        2      3           4          5           6        7           8   9       10      11  
predSubCleanedNoSpace = [initArea/1e6, temp, tempChange, precip, precipChange, glArea, glWidth/1e3, hterm, dbdz, b2010, cumMassBal ]; % remove lat, lon, elev, coastdist. Going to see how other variables cluster by these groups

predSubCleanedZeroMin = bsxfun(@minus,predSubCleaned,nanmin(predSubCleaned));
predSubCleanedRange1 = bsxfun(@rdivide,predSubCleanedZeroMin,nanmax(predSubCleanedZeroMin));

predSubCleanedZeroMean = bsxfun(@minus,predSubCleaned,nanmean(predSubCleaned));
predSubCleanedDivStd = bsxfun(@rdivide,predSubCleanedZeroMean,nanstd(predSubCleanedZeroMean));

predSubCleanedNoSpaceZeroMin = bsxfun(@minus,predSubCleanedNoSpace,nanmin(predSubCleanedNoSpace));
predSubCleanedNoSpaceRange1 = bsxfun(@rdivide,predSubCleanedNoSpaceZeroMin,nanmax(predSubCleanedNoSpaceZeroMin));

[coeff, score, trash, trash, varExp, trash] = pca(predSubCleanedRange1,'rows','complete');
[coeff2, score2, trash, trash, varExp2, trash] = pca(predSubCleanedDivStd,'rows','complete');
[coeff3, score3, trash, trash, varExp3, trash] = pca(predSubCleanedNoSpaceRange1,'rows','complete');


tableSortOrd = [6 7 8 9 5 10 11 12 13 14 15 1 2 3 4]; % index for sorting as data appear on table in paper
tableSortOrdNoSpace =  [2 3 4 5 6 7 8 9 10 11 1];

scoreWithLatLon = [lat,long,score];
%csvwrite('pca_loadings_withLatLon.csv',scoreWithLatLon)
%csvwrite('pca_loadings_divStd_sortForTable.csv',coeff2(tableSortOrd,:))
%csvwrite('pca_loadings_divStd_sortForTable.csv',coeff2(tableSortOrd,:))
csvwrite('pca_loadings_noSpace_sortForTable.csv',coeff3(tableSortOrdNoSpace,:))
%csvwrite('pca_varExp_divStd.csv',varExp2')

nums = [1:size(predSubCleaned,2)]';

[trash,sortVarInd] = sort(varExp,'descend');
[nums(sortVarInd),varExp(sortVarInd)]

pcNum = 1;

[trash,scoreSortInd] = sort(coeff(:,pcNum),'descend');
[trash,scoreSortInd2] = sort(coeff2(:,pcNum),'descend');

[nums(scoreSortInd), coeff(scoreSortInd,pcNum)]
[nums(scoreSortInd2), coeff(scoreSortInd2,pcNum)]
%                    1     2     3       4            5         6        7          8         9          10        11         12   13      14     15                   
% predSubCleaned = [lat, long, elev, initArea/1e6, coastDist, temp, tempChange, precip, precipChange, glArea, glWidth/1e3, hterm, dbdz, b2010, cumMassBal ];
varCell = {'lat','lon','elev','initArea','coastDist','temp','d temp','precip','d precip','gl area','gl width','h term','dbdz','b2010','bCum'};

% trying to interpret meanings of PCs. Just considering variables with >0.25 loading on that PC
% pc1 ~ just x/y location? = long (0.55), near term ice thick (-0.39), lat (-0.36), temp change (-.26)
% pc2 ~ continentality? = db/dz (.56), elev (-0.52), coast dist (-0.33), precip (0.32), temp (0.28)
% pc3 ~ glacier size = glacier area (.62), near terminus thickness (.53), longitude (.40), elev (.26)
% pc4 ~ climate / climate change = temp change (0.63), cum. mass bal (-0.36), 2010s mass balance (-0.35), precip (0.29), dbdz (-0.28), temp (-.28) 
% pc5 ~ also climate/change? = precip (0.62), 2010s mass bal (0.43), cumulative mass balance (0.37), temp (-0.29), elev (0.27)

predAndProglExistIndNoPCNan = predAndProglExistInd == 1 & isnan(score(:, 1)) == 0;
predAndDammedExistIndNoPCNan  = predAndDammedExistInd == 1 & isnan(score(:, 1)) == 0;

[pc1ProglAbsTau, pc1ProglAbsP] = corr(score(predAndProglExistIndNoPCNan, 1),dArea(predAndProglExistIndNoPCNan),'type','kendall');
[pc2ProglAbsTau, pc2ProglAbsP] = corr(score(predAndProglExistIndNoPCNan, 2),dArea(predAndProglExistIndNoPCNan),'type','kendall');
[pc3ProglAbsTau, pc3ProglAbsP] = corr(score(predAndProglExistIndNoPCNan, 3),dArea(predAndProglExistIndNoPCNan),'type','kendall');
[pc4ProglAbsTau, pc4ProglAbsP] = corr(score(predAndProglExistIndNoPCNan, 4),dArea(predAndProglExistIndNoPCNan),'type','kendall');

[pc1ProglRelTau, pc1ProglRelP] = corr(score(predAndProglExistIndNoPCNan, 1),dAreaRelStart(predAndProglExistIndNoPCNan),'type','kendall');
[pc2ProglRelTau, pc2ProglRelP] = corr(score(predAndProglExistIndNoPCNan, 2),dAreaRelStart(predAndProglExistIndNoPCNan),'type','kendall');
[pc3ProglRelTau, pc3ProglRelP] = corr(score(predAndProglExistIndNoPCNan, 3),dAreaRelStart(predAndProglExistIndNoPCNan),'type','kendall');
[pc4ProglRelTau, pc4ProglRelP] = corr(score(predAndProglExistIndNoPCNan, 4),dAreaRelStart(predAndProglExistIndNoPCNan),'type','kendall');

[pc1DamAbsTau, pc1DamAbsP] = corr(score(predAndDammedExistIndNoPCNan, 1),dArea(predAndDammedExistIndNoPCNan),'type','kendall');
[pc2DamAbsTau, pc2DamAbsP] = corr(score(predAndDammedExistIndNoPCNan, 2),dArea(predAndDammedExistIndNoPCNan),'type','kendall');
[pc3DamAbsTau, pc3DamlAbsP] = corr(score(predAndDammedExistIndNoPCNan, 3),dArea(predAndDammedExistIndNoPCNan),'type','kendall');
[pc4DamAbsTau, pc4DamAbsP] = corr(score(predAndDammedExistIndNoPCNan, 4),dArea(predAndDammedExistIndNoPCNan),'type','kendall');

[pc1DamRelTau, pc1DamRelP] = corr(score(predAndDammedExistIndNoPCNan, 1),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');
[pc2DamRelTau, pc2DamRelP] = corr(score(predAndDammedExistIndNoPCNan, 2),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');
[pc3DamRelTau, pc3DamlRelP] = corr(score(predAndDammedExistIndNoPCNan, 3),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');
[pc4DamRelTau, pc4DamRelP] = corr(score(predAndDammedExistIndNoPCNan, 4),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');


[mPC2AbsP, bPC2AbsP] = TheilSen([score(predAndProglExistIndNoPCNan,2),dArea(predAndProglExistIndNoPCNan)/1e6]);
[mPC2AbsD, bPC2AbsD] = TheilSen([score(predAndDammedExistIndNoPCNan,2),dArea(predAndDammedExistIndNoPCNan)/1e6]);

[mPC3AbsP, bPC3AbsP] = TheilSen([score(predAndProglExistIndNoPCNan,3),dArea(predAndProglExistIndNoPCNan)/1e6]);
[mPC3AbsD, bPC3AbsD] = TheilSen([score(predAndDammedExistIndNoPCNan,3),dArea(predAndDammedExistIndNoPCNan)/1e6]);

[mPC2RelP, bPC2RelP] = TheilSen([score(predAndProglExistIndNoPCNan,2),dAreaRelStart(predAndProglExistIndNoPCNan)*1e2]);
[mPC2RelD, bPC2RelD] = TheilSen([score(predAndDammedExistIndNoPCNan,2),dAreaRelStart(predAndDammedExistIndNoPCNan)*1e2]);



xl1 = [-.65,.6];
xl2 = [-.4,.8];

figure(293)
clf
orient landscape
subplot(2,1,1)
hold on
scatter(score(predAndProglExistIndNoPCNan, 2),dArea(predAndProglExistIndNoPCNan)/1e6,20,dbdz(predAndProglExistIndNoPCNan),'filled')
scatter(score(predAndDammedExistIndNoPCNan, 2),dArea(predAndDammedExistIndNoPCNan)/1e6,25,dbdz(predAndDammedExistIndNoPCNan),'linewidth',1.5,'marker','d')
plot(xl1,xl1.*mPC2AbsP + bPC2AbsP,'b','linewidth',2);
plot(xl1,xl1.*mPC2AbsD + bPC2AbsD,'r','linewidth',2);

subplot(2,1,2)
hold on
scatter(score(predAndProglExistIndNoPCNan, 2),dAreaRelStart(predAndProglExistIndNoPCNan)*1e2,20,dbdz(predAndProglExistIndNoPCNan),'filled')
scatter(score(predAndDammedExistIndNoPCNan, 2),dAreaRelStart(predAndDammedExistIndNoPCNan)*1e2,25,dbdz(predAndDammedExistIndNoPCNan),'linewidth',1.5,'marker','d')
plot(xl1,xl1.*mPC2RelP + bPC2RelP,'b','linewidth',2);
plot(xl1,xl1.*mPC2RelD + bPC2RelD,'r--','linewidth',1);

subplot(2,1,1)
ylim(absLim)
%caxis([0 200])
c = colorbar();
ylabel(c,'db/dz [m w.e. (100 m)^{-1}]')
ylabel('Absolute area change, \DeltaA [km^2]','fontsize',16)
grid on
plot([-.65,.6],[0 0],'k:')
xlim([-.65,.6])
%xlabel('PC2 score [-]')


subplot(2,1,2)
ylim(relLim)
%caxis([0 200])
c = colorbar();
ylabel(c,'db/dz [m w.e. (100 m)^{-1}]')
ylabel('Relative area change, \DeltaA/A_i [%]','fontsize',16)
grid on
plot([-.65,.6],[0 0],'k:')
xlim([-.65,.6])
xlabel('PC2 score [-]')

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('pc2_AbsAndRelChange.pdf','-dpdf','-r300')
%print('pc2_AbsAndRelChange.png','-dpng','-r300')


figure(294)

subplot(2,1,2)
hold on
scatter(score(predAndProglExistIndNoPCNan, 3),dArea(predAndProglExistIndNoPCNan)/1e6,20,glArea(predAndProglExistIndNoPCNan),'filled')
scatter(score(predAndDammedExistIndNoPCNan, 3),dArea(predAndDammedExistIndNoPCNan)/1e6,25,glArea(predAndDammedExistIndNoPCNan),'linewidth',1.5,'marker','d')

plot(xl2,xl2.*mPC3AbsP + bPC3AbsP,'b','linewidth',2);
plot(xl2,xl2.*mPC3AbsD + bPC3AbsD,'r--','linewidth',1);

subplot(2,1,2)
ylim(absLim)
%caxis([0 200])
c = colorbar();
ylabel(c,'Glacier area [km^2]')
ylabel('Area change [km^2]')
grid on
plot([-.4,.8],[0 0],'k:')
xlim([-.4,.8])
xlabel('PC3 score [-]')



% for the pca where data are demeaned and then div by stdev
[pc1ProglAbsTau2, pc1ProglAbsP2] = corr(score2(predAndProglExistIndNoPCNan, 1),dArea(predAndProglExistIndNoPCNan),'type','kendall');
[pc2ProglAbsTau2, pc2ProglAbsP2] = corr(score2(predAndProglExistIndNoPCNan, 2),dArea(predAndProglExistIndNoPCNan),'type','kendall');
[pc3ProglAbsTau2, pc3ProglAbsP2] = corr(score2(predAndProglExistIndNoPCNan, 3),dArea(predAndProglExistIndNoPCNan),'type','kendall');
[pc4ProglAbsTau2, pc4ProglAbsP2] = corr(score2(predAndProglExistIndNoPCNan, 4),dArea(predAndProglExistIndNoPCNan),'type','kendall');

[pc1ProglRelTau, pc1ProglRelP] = corr(score(predAndProglExistIndNoPCNan, 1),dAreaRelStart(predAndProglExistIndNoPCNan),'type','kendall');
[pc2ProglRelTau, pc2ProglRelP] = corr(score(predAndProglExistIndNoPCNan, 2),dAreaRelStart(predAndProglExistIndNoPCNan),'type','kendall');
[pc3ProglRelTau, pc3ProglRelP] = corr(score(predAndProglExistIndNoPCNan, 3),dAreaRelStart(predAndProglExistIndNoPCNan),'type','kendall');
[pc4ProglRelTau, pc4ProglRelP] = corr(score(predAndProglExistIndNoPCNan, 4),dAreaRelStart(predAndProglExistIndNoPCNan),'type','kendall');

[pc1DamAbsTau, pc1DamAbsP] = corr(score(predAndDammedExistIndNoPCNan, 1),dArea(predAndDammedExistIndNoPCNan),'type','kendall');
[pc2DamAbsTau, pc2DamAbsP] = corr(score(predAndDammedExistIndNoPCNan, 2),dArea(predAndDammedExistIndNoPCNan),'type','kendall');
[pc3DamAbsTau, pc3DamlAbsP] = corr(score(predAndDammedExistIndNoPCNan, 3),dArea(predAndDammedExistIndNoPCNan),'type','kendall');
[pc4DamAbsTau, pc4DamAbsP] = corr(score(predAndDammedExistIndNoPCNan, 4),dArea(predAndDammedExistIndNoPCNan),'type','kendall');

[pc1DamRelTau, pc1DamRelP] = corr(score(predAndDammedExistIndNoPCNan, 1),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');
[pc2DamRelTau, pc2DamRelP] = corr(score(predAndDammedExistIndNoPCNan, 2),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');
[pc3DamRelTau, pc3DamlRelP] = corr(score(predAndDammedExistIndNoPCNan, 3),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');
[pc4DamRelTau, pc4DamRelP] = corr(score(predAndDammedExistIndNoPCNan, 4),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');

[mPC2AbsP2, bPC2AbsP2] = TheilSen([score2(predAndProglExistIndNoPCNan,2),dArea(predAndProglExistIndNoPCNan)/1e6]);
[mPC2AbsD2, bPC2AbsD2] = TheilSen([score2(predAndDammedExistIndNoPCNan,2),dArea(predAndDammedExistIndNoPCNan)/1e6]);


% see both pca methods
figure(294)
clf
orient landscape
subplot(2,1,1)
hold on
scatter(score(predAndProglExistIndNoPCNan, 2),dArea(predAndProglExistIndNoPCNan)/1e6,20,'b','filled')
scatter(score2(predAndProglExistIndNoPCNan, 2),dArea(predAndProglExistIndNoPCNan)/1e6,20,'r','filled')
plot(xl1,xl1.*mPC2AbsP + bPC2AbsP,'b','linewidth',2);
plot(xl1,xl1.*mPC2AbsP2 + bPC2AbsP2,'r','linewidth',2);


%%
% see both pca methods
figure(310)
clf
orient landscape
subplot(2,1,1)
hold on
scatter(score(predAndProglExistIndNoPCNan, 1),score(predAndProglExistIndNoPCNan, 2),20,log10(dArea(predAndProglExistIndNoPCNan)),'filled')
xl = xlim;
yl = ylim;
xlim(xl)
ylim(yl)
plot(xl,[0 0],'k:')
plot([0 0],yl,'k:')
grid on
xlabel('PC1 score [-]','fontsize',16)
ylabel('PC2 score [-]','fontsize',16)
title('color = log10(dArea)')

subplot(2,1,2)
hold on
scatter(score(predAndProglExistIndNoPCNan, 2),score(predAndProglExistIndNoPCNan, 3),20,log10(dArea(predAndProglExistIndNoPCNan)),'filled')
xl = xlim;
yl = ylim;
plot(xl,[0 0],'k:')
plot([0 0],yl,'k:')
xlim(xl)
ylim(yl)
grid on
xlabel('PC2 score [-]','fontsize',16)
ylabel('PC3 score [-]','fontsize',16)
%colorbar()

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('pc1vs2_2vs3_withSpace_absAreaChange)withColorbar.pdf','-dpdf','-r300')


%%
figure(300)
clf
orient landscape
hold on
%scatter(score3(predAndProglExistIndNoPCNan,1),score3(predAndProglExistIndNoPCNan,2),coastDist(predAndProglExistIndNoPCNan)*1,hterm(predAndProglExistIndNoPCNan),'filled')
scatter(score3(predAndProglExistIndNoPCNan,1),score3(predAndProglExistIndNoPCNan,2),20,log10(dArea(predAndProglExistIndNoPCNan)),'filled')
%caxis([0 10])
colorbar()
%scatter(score3(predAndDammedExistIndNoPCNan,1),score3(predAndDammedExistIndNoPCNan,2),20,'rd','filled')

% for the pca without spatial data
[pc1ProglAbsTau3, pc1ProglAbsP3] = corr(score3(predAndProglExistIndNoPCNan, 1),dArea(predAndProglExistIndNoPCNan),'type','kendall');
[pc2ProglAbsTau3, pc2ProglAbsP3] = corr(score3(predAndProglExistIndNoPCNan, 2),dArea(predAndProglExistIndNoPCNan),'type','kendall');
[pc3ProglAbsTau3, pc3ProglAbsP3] = corr(score3(predAndProglExistIndNoPCNan, 3),dArea(predAndProglExistIndNoPCNan),'type','kendall');
[pc4ProglAbsTau3, pc4ProglAbsP3] = corr(score3(predAndProglExistIndNoPCNan, 4),dArea(predAndProglExistIndNoPCNan),'type','kendall');


[pc1DamAbsTau3, pc1DamAbsP3] = corr(score3(predAndDammedExistIndNoPCNan, 1),dArea(predAndDammedExistIndNoPCNan),'type','kendall');
[pc2DamAbsTau3, pc2DamAbsP3] = corr(score3(predAndDammedExistIndNoPCNan, 2),dArea(predAndDammedExistIndNoPCNan),'type','kendall');
[pc3DamAbsTau3, pc3DamlAbsP3] = corr(score3(predAndDammedExistIndNoPCNan, 3),dArea(predAndDammedExistIndNoPCNan),'type','kendall');
[pc4DamAbsTau3, pc4DamAbsP3] = corr(score3(predAndDammedExistIndNoPCNan, 4),dArea(predAndDammedExistIndNoPCNan),'type','kendall');


[pc1DamRelTau3, pc1DamRelP3] = corr(score3(predAndDammedExistIndNoPCNan, 1),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');
[pc2DamRelTau3, pc2DamRelP3] = corr(score3(predAndDammedExistIndNoPCNan, 2),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');
[pc3DamRelTau3, pc3DamlRelP3] = corr(score3(predAndDammedExistIndNoPCNan, 3),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');
[pc4DamRelTau3, pc4DamRelP3] = corr(score3(predAndDammedExistIndNoPCNan, 4),dAreaRelStart(predAndDammedExistIndNoPCNan),'type','kendall');


[mPC2AbsP3, bPC2AbsP3] = TheilSen([score3(predAndProglExistIndNoPCNan,2),dArea(predAndProglExistIndNoPCNan)/1e6]);
[mPC2AbsD3, bPC2AbsD3] = TheilSen([score3(predAndDammedExistIndNoPCNan,2),dArea(predAndDammedExistIndNoPCNan)/1e6]);

%% testing withheld variable correlations on pc axes
[coastDistPc1Tau,coastDistPc1P] = corr(coastDist(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 1),'type','kendall');
[coastDistPc2Tau,coastDistPc2P] = corr(coastDist(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 2),'type','kendall');
[coastDistPc3Tau,coastDistPc3P] = corr(coastDist(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 3),'type','kendall');

[elevPc1Tau,elevPc1P] = corr(elev(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 1),'type','kendall');
[elevPc2Tau,elevPc2P] = corr(elev(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 2),'type','kendall');
[elevPc3Tau,elevPc3P] = corr(elev(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 3),'type','kendall');

[latPc1Tau,latPc1P] = corr(lat(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 1),'type','kendall');
[latPc2Tau,latPc2P] = corr(lat(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 2),'type','kendall');
[latPc3Tau,latPc3P] = corr(lat(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 3),'type','kendall');

[longPc1Tau,longPc1P] = corr(long(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 1),'type','kendall');
[longPc2Tau,longPc2P] = corr(long(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 2),'type','kendall');
[longPc3Tau,longPc3P] = corr(long(predAndProglExistIndNoPCNan),score3(predAndProglExistIndNoPCNan, 3),'type','kendall');


%%
xl3 = [-0.6 0.8];

figure(301)
clf
orient landscape
hold on
scatter(score3(predAndProglExistIndNoPCNan, 2),dArea(predAndProglExistIndNoPCNan)/1e6,25,coastDist(predAndProglExistIndNoPCNan),'filled')
scatter(score3(predAndDammedExistIndNoPCNan, 2),dArea(predAndDammedExistIndNoPCNan)/1e6,20,coastDist(predAndDammedExistIndNoPCNan),'marker','d','linewidth',1.5)
plot(xl3,xl3.*mPC2AbsP3 + bPC2AbsP3,'b','linewidth',2);
plot(xl3,xl3.*mPC2AbsD3 + bPC2AbsD3,'r','linewidth',1);
plot(xl3,[0 0],'k:','linewidth',1);
plot([0 0],absLim,'k:','linewidth',1);

caxis([0 300])
cb = colorbar();
ylabel(cb,'Distance from coast [km]','fontsize',16)
grid on
set(gca,'fontsize',12)
ylabel('Absolute area change, \DeltaA [km^2]','fontsize',16)
xlabel('PC2 score [-]','fontsize',16)
ylim(absLim)
xlim(xl3)

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
%print('pc2_noSpace_AbsAndRelChange.pdf','-dpdf','-r300')


%% run cross-correlation for input variables

if 0
    corrMatxTau = nan(length(nums),length(nums));
    corrMatxP = nan(length(nums),length(nums));


    figure(295)
    clf
    orient landscape

    for i = 2:length(nums)
        varNow =predSubCleaned(:,i);

        firstDataInd = ~isnan(varNow);

        for j = 1:length(nums)
           otherVarNow = predSubCleaned(:,j);
           secondDataInd = ~isnan(otherVarNow);
           bothDataInd = firstDataInd == 1 & secondDataInd == 1;

           [tauNow, pNow] = corr(varNow(bothDataInd),otherVarNow(bothDataInd),'type','kendall');
           corrMatxTau(i,j) = tauNow;
           corrMatxP(i,j) = pNow;

           [mNow, bNow] = TheilSen([varNow(bothDataInd),otherVarNow(bothDataInd)]);

           %subplot(15,15,(i*15)+j)


           figure(295)
           clf
           orient landscape
           hold on

           scatter(varNow(bothDataInd)',otherVarNow(bothDataInd)','filled')
           if pNow < 0.05
            plot(xlNow,xlNow.*mNow + bNow,'b','linewidth',2);
           end
            xlNow = xlim();
            ylNow = ylim();
            text(xlNow(1),ylNow(2),['tau = ',num2str(tauNow)])       
           text(xlNow(2),ylNow(2),['p = ',num2str(pNow)])               
           title(['x = ',cell2mat(varCell(i)),'; y = ',cell2mat(varCell(j))])

           %pause()

           h = gcf;
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);
            outFn = ['crossCorrPlots_', cell2mat(varCell(i)), ' _vs_', cell2mat(varCell(j)) '.pdf'];
            print(outFn,'-dpdf','-r300')

        end

    end

end

%csvwrite('controlVariables_crossCorrelation_tau.csv',corrMatxTau)
%csvwrite('controlVariables_crossCorrelation_p.csv',corrMatxP)


%% trying to figure out function to automatically point to off plot data, didn't finish
if 0
xData = dbdz(predAndProglExistInd);
yData = dAreaRelStart(predAndProglExistInd)*100;
xLims = xl2;
yLims = relLim;
yRange = yLims(2) - yLims(1);
yOffset = yRange/40;


offYind = yData < yLims(1) | yData > yLims(2);
xForOut = xData(offYind);
yOff = yData(offYind);

belowYind = yOff < yLims(1);
aboveYind = yOff > yLims(2);

yForOut = nan(length(yOff),1);
yForOut(aboveYind) = yLims(2) - yOffset;
yForOut(belowYind) = yLims(1) + yOffset;

figure(99)
clf
hold on
scatter(xData,yData)
scatter(xData(outYind),yData(outYind),'r','filled')
anArrow = annotation('arrow','color','b');
anArrow.Parent = gca;
anArrow.Position = [xForOut(1),yLims(2)-yRange/40,0,yRange/40];

ylim(yLims)

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print('test.pdf','-dpdf','-r300')
print('test.png','-dpng','-r300')

end