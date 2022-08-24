
% remake figures from Reza for WE43-T6 paper part 2, matching the style of
% my other figures

working_dir = 'C:\Users\chenz\Downloads\FiguresOrigin';

figure;
fig_position = get(gcf,'position');
axis_position = get(gca,'position');
close;

%% fig 16a
close all;
fig_name = 'fig 16a';
open(fullfile(working_dir,[fig_name,'.fig']));

set(gcf,'WindowState','normal');
set(gcf,'position',fig_position);
set(gca,'fontsize',16,'ylim',[0 1.4],'YTick',[0:0.2:1.2],'LineWidth',1);
legend('fontsize',16);
print(fullfile(working_dir, [fig_name,'.tiff']), '-dtiff');

%% fig 16b
close all;
fig_name = 'fig 16b';
open(fullfile(working_dir,[fig_name,'.fig']));

set(gcf,'WindowState','normal');
set(gcf,'position',fig_position);
set(gca,'fontsize',16,'ylim',[0 1.4],'YTick',[0:0.2:1.2],'LineWidth',1);
legend('fontsize',16);
print(fullfile(working_dir, [fig_name,'.tiff']), '-dtiff');

%% fig 17a
close all;
fig_name = 'fig 17a';
open(fullfile(working_dir,[fig_name,'.fig']));

set(gcf,'WindowState','normal');
set(gcf,'position',fig_position);
set(gca,'fontsize',16,'LineWidth',1,'YTick',[0:0.1:1]);
hold on; plot(nan,nan,'.w');
legend({'SEM-DIC Experiment'},'fontsize',14);
print(fullfile(working_dir, [fig_name,'.tiff']), '-dtiff','-r300');

%% fig 17b
close all;
fig_name = 'fig 17b';
open(fullfile(working_dir,[fig_name,'.fig']));

set(gcf,'WindowState','normal');
set(gcf,'position',fig_position);
set(gca,'fontsize',16,'LineWidth',1,'YTick',[0:0.1:1]);
hold on; plot(nan,nan,'.w');
legend('PRISMS-Plasticity Simulation','fontsize',14);
print(fullfile(working_dir, [fig_name,'.tiff']), '-dtiff','-r300');

%% fig 18
close all;
fig_name = 'fig 18';
open(fullfile(working_dir,[fig_name,'.fig']));

set(gcf,'WindowState','normal');
set(gcf,'position',fig_position);
set(gca,'fontsize',16,'LineWidth',1,'YTick',[0:0.1:1]);
legend('fontsize',16);
print(fullfile(working_dir, [fig_name,'.tiff']), '-dtiff');

%% fig 19a
close all;
fig_name = 'fig 19a';
open(fullfile(working_dir,[fig_name,'.fig']));

set(gcf,'WindowState','normal');
set(gcf,'position',fig_position);
set(gca,'fontsize',16,'LineWidth',1,'YTick',[0:0.1:1]);
hold on; plot(nan,nan,'.w');
legend({'SEM-DIC Experiment'},'fontsize',14);
print(fullfile(working_dir, [fig_name,'.tiff']), '-dtiff');

%% fig 19b
close all;
fig_name = 'fig 19b';
open(fullfile(working_dir,[fig_name,'.fig']));

set(gcf,'WindowState','normal');
set(gcf,'position',fig_position);
set(gca,'fontsize',16,'LineWidth',1,'YTick',[0:0.1:1]);
hold on; plot(nan,nan,'.w');
legend('PRISMS-Plasticity Simulation','fontsize',14);
print(fullfile(working_dir, [fig_name,'.tiff']), '-dtiff');

%%
close all;


