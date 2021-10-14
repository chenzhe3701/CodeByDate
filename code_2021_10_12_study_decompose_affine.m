% How do we estimate global strain? epsilon from deformation gradient, or
% decompose transformation into steps of transformation?
% Seems that the results are similar if deformation is small.
% But epsilon from deformation gradient is more reasonable.
% 2021-10-13

% example, study the relationship between deformation gradient and decompose affine 2d  
close all;
clc;

%% example from Mg4Al_U2, iE=0 and iE=2
cpFrom = [
    97.7540   71.1240;
  474.8690  107.3170;
   88.9200  473.3950];
cpTo = [
      101.2160   67.8020;
  470.4430  103.0320;
   95.3930  459.6990];

%% concrete example of randomly generated cases
cpFrom = [
    1.9643   -4.0765;
   16.2533   -4.9218;
   -3.6985   19.2311];
cpTo = [
    1.5557   -3.9118;
   22.2292    1.3177;
    0.3121   16.2650];

%% randomly generated example, more extreme cases
cpFrom = [0 0; 20,0; 0,20] + (rand(3,2)-0.5)*10;
cpTo = [0 0; 20,0; 0,20] + (rand(3,2)-0.5)*10;

%
tform = fitgeotrans(cpFrom, cpTo, 'affine');

ptFrom = [cpFrom(1,:),1]';
ptTo = [cpTo(1,:),1]';

M = tform.T';
M * ptFrom;
ptTo;
% we should have: M * ptFrom = ptTo, i.e., x = MX

[T,R,Z,S,A,epsilon] = decompose_affine2d(M);
Tm = [1,0,T(1); 0,1,T(2); 0,0,1];
Rm = R;
Zm = [Z(1), 0, 0; 0, Z(2), 0; 0, 0, Z(3)];
Sm = [1, S(1), S(2); 0, 1, S(3); 0 0 1];

RZS = Rm * Zm * Sm;
A;
% we should have Rm * Zm * Sm = A(1:3,1:3)

% If similar to 2D-DIC, we want to have an F=dx/dX=du/dX+I that has the relationship: x = F*X  
% according to the definition of deformation gradient:
% dx/dX = M(1,1), dx/dY = M(1,2)
% dy/dX = M(2,1), dy/dY = M(2,2)
% F = M(1:2,1:2)
F = M(1:2,1:2);
epsilon = (F'*F-eye(2))/2
Zm-1

% plot step by step
close all;
X = [cpFrom(1:3,:)'; ones(1,3)]; % convert to column format
x = [cpTo(1:3,:)'; ones(1,3)];

X_1 = Sm * X;
X_2 = Zm * X_1;
X_3 = Rm * X_2;
X_4 = Tm * X_3;

colors = parula(5);
colors(5,:) = [1 0 0];
figure; hold on;
plot(X(1,[1:3,1]), X(2,[1:3,1]), '-o', 'color','k', 'linewidth',2);
plot(X_1(1,[1:3,1]), X_1(2,[1:3,1]), '-o', 'color',colors(2,:));
plot(X_2(1,[1:3,1]), X_2(2,[1:3,1]), '-o', 'color',colors(3,:));
plot(X_3(1,[1:3,1]), X_3(2,[1:3,1]), '-o', 'color',colors(4,:));
plot(X_4(1,[1:3,1]), X_4(2,[1:3,1]), '-o', 'color',colors(5,:), 'linewidth',2);
plot(x(1,[1:3,1]), x(2,[1:3,1]), '-o', 'color','r');
str = {'Reference', ['Shear xy=',num2str(S(1),3)], ['Zoom x=',num2str(Z(1),3),' y=',num2str(Z(2),3)], 'Rotate', 'Translate'};
set(gca,'fontsize',14);
xlabel('X');
ylabel('Y');
% add_legend(str);
legend(str,'location','best');



%% Which estimate is better? Zm vs. epsilon?