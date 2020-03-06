% code 2020_03_01
% summary the relationship between HAuCl4, Na3Citrate, and Au diameter 
xy = [1, 16;
    0.75, 25;
    0.5, 41;
    0.3, 72;
    0.21, 98;
    0.67, 40;
    0.43, 60;
    1, 25;
    0.63, 60;
    0.83, 45;
    0.44, 55;];
x = xy(:,1);
y = xy(:,2);

figure; hold on;
plot(x,y,'o');
plot(x(1:5),y(1:5),'or','markerfacecolor','r')
set(gca,'yScale','log')