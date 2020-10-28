% study how to hedge for the game
close all;

invest_A = 1000;
invest_B = 1500:10:2500;

odd_A = 0.4;
odd_B = 1.7;

return_A = invest_A ./ odd_A * 0.9 - invest_B;   % if A wins
return_B = invest_B ./ odd_B * 0.9 - invest_A;  % if B wins

figure;
hold on;
plot(invest_B, return_A);
plot(invest_B, return_B);
plot(invest_B, zeros(size(invest_B)));
legend('A wins', 'B wins');


