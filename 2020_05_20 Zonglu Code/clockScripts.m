figure;
for hour = 0:12
    for minute = 0:59
        plotTime(hour, minute);
        pause(0.05);
        clf;
    end
end