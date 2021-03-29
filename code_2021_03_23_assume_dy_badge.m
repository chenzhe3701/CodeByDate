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