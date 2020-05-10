function [] = plotTime(hour, minute)

t = 0 : 2*pi/100 : 2*pi;
x = cos(t);
y = sin(t);

hold on;
plot(x,y,'b')
axis equal;

w = 0.05;
L = 0.9;


minuteHand = [-w/2 w/2 w/2 -w/2 -w/2; 0 0 L L 0];

g = [cosd(-minute/60*360), -sind(-minute/60*360);
    sind(-minute/60*360), cosd(-minute/60*360)];
rotatedMinuteHand = g * minuteHand;

plot(rotatedMinuteHand(1,:), rotatedMinuteHand(2,:),'g');


hourHand = [-w/2 w/2 w/2 -w/2 -w/2; 0 0 L/2 L/2 0];

g = [cosd(-hour/12*360 - minute/60*360/12), -sind(-hour/12*360 - minute/60*360/12);
    sind(-hour/12*360 - minute/60*360/12), cosd(-hour/12*360 - minute/60*360/12)];
rotatedHourHand = g * hourHand;

plot(rotatedHourHand(1,:), rotatedHourHand(2,:), 'r');

end