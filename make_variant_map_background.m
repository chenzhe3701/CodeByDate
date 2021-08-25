function [] = make_variant_map_background()

colors = parula(7);
colors(1,:) = [0.25, 0.25, 0.25];
colormap(colors);

a = findall(gcf,'type','Axes');
set(a,'xTick',[], 'yTick',[], 'fontsize',16);
title('');

caxis([-0.5, 6.5]);
c = findall(gcf,'type','Colorbar');
set(c, 'limits', [0.5,6.5]);

