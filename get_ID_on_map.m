
% helpler function to get grain ID from a map.
% 2021-01-23

function id = get_ID_on_map(ID)

h = drawpoint;
customWait(h);
pos = round(h.Position);

indr = pos(2);
indc = pos(1);

id = ID(indr, indc);

display(['grain ID = ',num2str(id)]);

