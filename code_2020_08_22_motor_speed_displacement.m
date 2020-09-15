speed = [0.5,1,2,4,8];

for ii = 1:length(speed)
  fid = fopen(['data\speed_displacement_measure\measurement_speed_',num2str(speed(ii)),'.lvm']);
  t1 = textscan(fid,'%s', 5, 'Delimiter', '\t');
  t2 = textscan(fid, '%f %f %f %f %f');
  fclose(fid); 
  data = cell2mat(t2);
  data(1,:) = [];
  data(end,:) = [];
  
  time(ii) = data(end,1) - data(1,1);
  displacement(ii) = abs(data(end,2) - data(1,2));
  md = fitlm(data(:,1),data(:,2));
  rate(ii) = abs(md.Coefficients.Estimate(2));
  voltage(ii) = speed(ii)/10*5;
  
end

md = fitlm(rate*1000,voltage);
volts_meter_micron_per_second = md.Coefficients.Estimate(2)    % volts per (um/s), 0.1571V for 1um/s

displacement_SEM = [1,2,2,2,2];
displacement_SEM_over_meter = mean(displacement_SEM./displacement);
volts_SEM_micro_per_second = volts_meter_micron_per_second * displacement_SEM_over_meter    % if read by SEM, 0.1542V for 1um/s