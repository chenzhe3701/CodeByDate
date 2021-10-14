% add legend, but keep figure height/weight the same
% 2021-10-12
function [] = add_legend(str)

set(gca,'unit','pixels');
pa_0 = get(gca,'position');
set(gca,'unit','normalized');

legend(str, 'location','eastoutside');

% if vertical position of axes does not change, it means we need to check for the horizontal (width)  
while(true)
    pf = get(gcf,'position');
    pf = pf + [0 0 1 0];
    set(gcf, 'position',pf);
    
    set(gca,'unit','pixels');
    pa = get(gca,'position');
    set(gca,'unit','normalized');
    if pa(2)>pa_0(2) | pa(3)>pa_0(3)
        break;
    end
end

end