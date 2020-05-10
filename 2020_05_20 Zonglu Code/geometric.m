function [xt, yt] = geometric(x, y, option)

switch option
    case 0
        xt = x;
        yt = y;
    case 1
        xt = x + 2;
        yt = y - 3;
    case 2
        xt = x * 0.5;
        yt = y * 2;
    case 3
        g = [cosd(-45), -sind(-45);
            sind(-45), cosd(-45)];
        temp = g * [x;y];
        xt = temp(1,:);
        yt = temp(2,:);
    case 4
        xt = y;
        yt = x;
end

end