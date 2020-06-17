% T1 column: {'iE','ID','e1','e5','e10','e50','e90','e95','e99','emean','basal_SF','pm_SF','py_SF','pyII_SF','twin_SF','activeTwin_SF','gDia','twinnedTF','tAF'};
% T2 column: {'iE','ID','iTwin','variant','variant_SF','vActiveTF','vPct','v1_activeTF','basal_SF','twin_SF'};

load('WE43_T6_C1_tableSummary_twinStats_selfEffect.mat','T','T2');

%%
disp(['========twin variant activated/not vs SF==============================']);
colororder('default');
C = colororder;
C = [.5 .5 .5; 1 0 0];
edges = -0.5:0.05:0.5;
for iE=2:5
    ind = (T2.iE==iE)&(T2.vActiveTF == 1);
    ind2 = (T2.iE==iE)&(T2.vActiveTF == 0);
    [N_t,~] = histcounts(T2.variant_SF(ind), edges);
    [N_nt,~] = histcounts(T2.variant_SF(ind2), edges);
    
    figure; hold on; disableDefaultInteractivity(gca);
    hbar = bar(edges(1:end-1)+0.025, [N_nt(:), N_t(:)], 1, 'stacked');
    set(gca,'fontsize',12, 'XTick',-0.5:0.1:0.5,'ylim',[0 1600]);
    xlabel('Twin Variant Schmid Factor');
    ylabel('Counts');
    
    yyaxis right;
    set(gca, 'ycolor', 'k','fontsize',16);
    ylabel('Percent');
    plot(-0.475:0.05:0.475, N_t./(N_t+N_nt),'-ko','linewidth',1.5);

    title(['iE = ',num2str(iE)],'fontweight','normal');
    legend({'Variants not twinned', 'Variants twinned','Percent of variants twinned'},'Location','northwest');
    
    hbar(1).FaceColor = [0 0 1];
    hbar(2).FaceColor = [1 0 0];
    
    % [for Plot used for paper, Fig 4a]:  X = Twin Variant Schmid Factor,  Y = Counts,  YY = Percent,  for All Twin Variants
    if iE==5
        title('All Twin Variants, \epsilon^G = -0.039');
        disp('row1: # twinned, row2: # not twinned, row3: pct');
        tt = [N_t; N_nt; N_t./(N_t+N_nt)];
        disp(array2table(tt));
    end
end
