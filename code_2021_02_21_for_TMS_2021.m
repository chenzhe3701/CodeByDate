% This can replace code for plotting, to use different color for the plots
% to show the effect of m' factor.
%
% To run the first part, look at code:
% script_paper_2020_part_2_analyze_twin_activity_summary

working_dir = 'E:\zhec umich Drive\0_temp_output';
close all;
pmf_all_exp = (N1+N2)./sum([N1(:);N2(:)]);  % probability mass function
pmf_basal_exp = N1./sum(N1(:));  % probability mass function, for active variants interacting with neighbor's basal
pmf_twin_exp = N2./sum(N2(:));  % probability mass function, for active variants interacting with neighbor's twin


% [Plot for paper, Fig 12a ]
figure; disableDefaultInteractivity(gca);
bar(xpos, [N1(:)+N2(:)], 1, 'stacked', 'facecolor', [0 0.4470 0.7410], 'facealpha', 1);
set(gca,'fontsize',16,'ylim',[0 170]);
xlabel('Representative m'' factor');
ylabel('Counts');
title('Active Variants','fontweight','normal');
print(fullfile(working_dir,'m'' active variant'), '-dtiff');

% [Plot for paper, Fig 12b ]
% (7) For: all the possible variants, look at their m' wrt all neighbor grains' MOST POSSIBLE slip or twin
N5 = histcounts(TT.mPrime, edges);
figure; disableDefaultInteractivity(gca);
bar(xpos, [N5(:)], 1, 'stacked','facecolor', 'k', 'facealpha', 0.8);
set(gca,'fontsize',16,'ylim',[0 7200]);
xlabel('Representative m'' factor');
ylabel('Counts');
title('All Possible Variants','fontweight','normal');
pmf_all = N5./sum(N5);
print(fullfile(working_dir,'m'' all variant'), '-dtiff');


% [Plot for paper, Fig 12c ]
% (8) relative strength, slip/twin in neighbor
figure; disableDefaultInteractivity(gca); hold on;
plot(xpos, pmf_all_exp, '-ko','linewidth',1.5,'MarkerSize',8, 'color', [0 0.4470 0.7410] );
plot(xpos, pmf_all, '--kd','linewidth',1.5,'MarkerSize',8);
set(gca,'ylim',[0 0.2]);
ylabel('Probability');
yyaxis right;
plot(xpos, pmf_all_exp./pmf_all, '-ro','linewidth',1.5,'MarkerSize',8);
set(gca,'ycolor','r', 'ylim',get(gca,'ylim').*[0,1]);
ylabel('Ratio');
xlabel('Representative m'' factor');
set(gca,'fontsize',16);
legend({'Probability, active variants','Probability, all possible variants','Ratio, active/all pissible variants'},'location','northwest','fontsize',17);
title('');
print(fullfile(working_dir,'m'' ratio'), '-dtiff');

%%