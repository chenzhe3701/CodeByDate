function [] = myfun()

l = 14;
w = 10;
h = 3.5;

L = [l, 0, 0];
W = [0, w, 0];
H = [0, 0, h];

pt{1} = [0,0,0];
pt{2} = pt{1} + W;
pt{3} = pt{2} + H;
pt{4} = pt{1} + H;

pt{5} = pt{1} + L;
pt{6} = pt{2} + L;
pt{7} = pt{3} + L;
pt{8} = pt{4} + L;

pt{9} = [8.5,0,0];
pt{10} = [5.5,w,0];
pt{11} = [6.5,w,h];
pt{12} = [9.5,0,h];

p1 = [5, w/2, 0];
pA = [7, w/2, 0];
p2 = [9, w/2, 0];

o = [-2.5, 0, -1];
ox = o + [2, 0, 0];
oy = o + [0, 4, 0];
oz = o + [0, 0, 2];

% right part of crystal has displacement s w.r.t left part of crystal
ra = 0.08;
rb = 0.16;
s = ra * (pt{9} - pt{10}) + rb * (pt{12} - pt{9});  % make sure the displacement is on this plane

pt{13} = pt{9} + s;
pt{14} = pt{10} + s;
pt{15} = pt{11} + s;
pt{16} = pt{12} + s;

pt{5} = pt{5} + s;
pt{6} = pt{6} + s;
pt{7} = pt{7} + s;
pt{8} = pt{8} + s;

pB = pA + s;
p2 = p2 + s;

% plot
close all;
figure('Position',[100,100,1200,600]);
view(22,20);
hold on;

plot_face({pt{1},pt{2},pt{10},pt{9}},'w', 0);
plot_face({pt{2},pt{3},pt{11},pt{10}},'w', 0);
plot_face({pt{5},pt{6},pt{7},pt{8}},'w', 0);
plot_face({pt{13},pt{14},pt{6},pt{5}},'w', 0);
plot_face({pt{14},pt{15},pt{7},pt{6}},'w', 0);

plot_face({pt{9},pt{10},pt{11},pt{12}},'k', 0.3);

plot_line(pt{1},pt{4},'--k');
plot_line(pt{4},pt{3},'--k');
plot_line(pt{4},pt{12},'--k');
plot_line(pt{15},pt{16},'--k');
plot_line(pt{13},pt{16},'--k');
plot_line(pt{16},pt{8},'--k');

plot_line(o,ox,'-k');
plot_line(o,oy,'-k');
plot_line(o,oz,'-k');
text(ox(1) + l/50, ox(2)+w/20, ox(3), 'x','color','k','fontsize',32);
text(oy(1), oy(2)+w/10, oy(3), 'y','color','k','fontsize',32);
text(oz(1), oz(2), oz(3)+h/10, 'z','color','k','fontsize',32);

plot_line(p1,pA,'-r');
plot_line(pA,pB,'-r');
plot_line(pB,p2,'-r');
text(p1(1)-l/15, p1(2)+w/15, p1(3), 'P_1','color','r','fontsize',32);
text(p2(1)+l/30, p2(2)+w/20, p2(3), 'P_2','color','r','fontsize',32);

plot3(p1(1), p1(2), p1(3), 'r.', 'markersize', 36, 'linewidth', 4);
plot3(p2(1), p2(2), p2(3), 'r.', 'markersize', 36, 'linewidth', 4);


set(gca,'ydir','reverse','zdir','reverse');
axis equal;
axis off;

end

function [] = plot_line(pt1,pt2,marker)

plot3([pt1(1),pt2(1)], [pt1(2),pt2(2)], [pt1(3),pt2(3)], marker, 'linewidth', 2);

end

function [] = plot_face(ptCell, color, facealpha)

m = [];
for ii = 1:length(ptCell)
    m = [m; ptCell{ii}];
end

fill3(m(:,1), m(:,2), m(:,3), color, 'facealpha', facealpha, 'LineWidth', 2);

end