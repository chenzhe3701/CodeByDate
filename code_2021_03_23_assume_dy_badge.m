% assumption for douyu badge
a(1,1) = 30;
a(1,2) = 100000;

exp_step = 10000;
experience = 100000;
for iL = 31:50
   level = iL;
   experience = experience + exp_step;
   exp_step = exp_step + 2000;
   a(iL-29,1) = iL;
   a(iL-29,2) = experience;
end
open a;

%% 2021-06-13 adjust rate
rate = [1.5; 1.5; 3.35; 4.9; 6.9; 11.9; 19.85; 24.8; 49.5; 1.5; 3.35]/100;
increase = [0.15; 0.14; 0.13; 0.12; 0.11; 0.1; 0.08; 0.06; 0.05; 0.2; 0.2];
rate_adjusted = rate .* (1+increase);

pay = [30; 20; 17; 15; 14; 12; 10; 5; 1; 30; 17];
unit_return = [2000; 1314; 500; 300; 200; 100; 50; 20; 2; 2000; 500];

expectation = rate .* unit_return ./ pay;
expectation_adjusted = rate_adjusted .* unit_return ./ pay;