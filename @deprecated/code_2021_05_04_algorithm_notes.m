% Need to check:
% (1) if we do segmentation of the twin map, will it get the same number of
% grains as in the original data?

% pixel level, at each iE, there are:
% (p1) current twin (directly from twin map)
% (p2) ever twinned up to this iE (obtain by accummulating twin map) 
% based on these two basic maps, we can get
% (p3, yellow) new twinned = currently twinned - ever twinned before = p1 - p2(iE-1);  
% (p4, blue) re-twinned = currently twinned & ever twinned bofore = p1 & p2(iE-1); 
% (p5, gray) detwinned = ever twinned before - currently twinned = p2(iE-1) - p1;
% p2 = p3 + p4 + p5
% (p6, black) completely detwinned area, a subset of p5, this should use
% grain level data
% 
% So, we have these pixel level maps:
% [1] p1 map, current twin map;
% [2] p2 map, ever twinned map; (p1 is a subset of p2)
% [3] p3 + p4 + p5 map, new/old/de-twin map; (p1=p3+p4, p5=p2-p1 <-> p2 = p3+p4+p5)
% 
%
% 
% grain level, we have the [1] ever twinned grain, [2] currently twinned grain   
%
% at each iE, the [current twin grains] can be categorized into:
% (g1) first time new twin, which contains only p3/new-twin pixel, no p4/re-twin, (no p5/de-twin)  
% (g2) completely detwin then new twinned (i.e., overlap with completely detwinned g4 at previous strain level, iE-1)
% (g3) mixed, which contains maybe p3/new-twin and at least p4/re-twin, (no p5/de-twin) 
%
% at each iE, we can GENERATE the map of the [ever twinned grains]:  
% (1) if a current twin grain does not overlap with any 'ever twin grain', 
% then add this grain to the 'ever-twinned map';
% (2) if a current twin grain overlap with an 'ever twin grain', merge it with
% the 'ever twin grain'
% So, ever twin grain = mixed + completely-detwinned
% (g4) mixed, contains any of p3/new-twin, p4/re-twin (p1/current-twin = p3+p4), p5/de-twin pixels
% (g5) completely detwinned, contains ONLY p5/detwinned
% ==> change g2 to 'ever completely detwin'? then we can check, if an ever
% completely detwinned grain currently have twin pixel, it is 'completely
% detwinned then retwinned grain'.


% ===> maybe
% current twin grain: 
% (g1) completely new twin:
%   Does not overlap with previous strain level twin(iE-1). 
%   Contains only new twin pixel at this strain level (iE)  
% (g3) evoling twin: overlap with prevoius strain level twin(iE-1)
% (g2) detwin then retwin twin: does not overlap with previous strain level
% twin, 



% grain wise:
% cycle-1-compression, iE=0:3
% (g1) new twin (contains only p3/new-twin pixel, no p4/re-twin, no p5/de-twin: yellow), 
% (g2) mix = growing twin (contains p3/new-twin + p4/re-twin, no p5/de-twin: yellow + blue) 
% cycle-1-tension, iE = 4:7
% (g1) new twin, should be very rare (only p3/new-twin pixel: yellow) 
% (g4) completely detwinned (only p5/detwin pixel: black) 
% (g2) currently twin part of a mixed/partially detwin grain (p4, maybe some p3,: blue, maybe some yellow)
% (g3) ever twin part of a mixed/partially detwin grain (p4+p5, maybe some p3: blue + gray, maybe some yellow)

 
%
% cycle 2, compression, iE=8:10
% (1) completely new twin (only new_twin pixel, p3, yellow)
% (2) completely detwinned, not retwinned (only detwin pixel, p6, black)
% (3) completely detwinned, then retwinned (detwin + retwin, p6+p4, black+blue)
% (4) mixed (detwinned + retwinned + new_twined, p5+p4+p3, gray+blue+yellow)
% cycle 2, tension, iE=11:13
% (1) completely detwinned (only detwinned pixel, p6, black) 
% (2) mixed (residual/re-twinned + detwinned pixel, p4+p5, blue + gray)
% (3) new twin grain, should be very rare (only new_twin pixel, p3, yellow)





% cycle-1-compression:
% (g1) new twin, which contains only p3/new-twin pixel, no p4/re-twin, no p5/de-twin
% (g2) mixed = growing twin, which contains both p3/new-twin and p4/re-twin, no p5/de-twin 
% cycle-2-tension:
% (1) mixed/partially detwinned grain, contains p4, maybe p3 pixels
% completely detwinned grain
% (2) new twin, contains only p3 pixel, and no p4, no p5
% cycle-2-compression:
% (1) new twin, contains only p3 pixel, and no p4, no p5; 
% (2) detwinne
