% recorded data point during insitu test for Mg4Al s1
close all;

um = [50, 100, 168, 284, 363, 442]
ue = [900 1800 7500 14900, 20000, 25000]
figure;
plot(um,ue,'-o')
set(gca,'xlim',[0,500],'ylim',[0 25000])
xlabel('disp, um');
ylabel('strain, ue');

md = fitlm(um(3:end),ue(3:end));
slope = md.Coefficients.Estimate(2)
