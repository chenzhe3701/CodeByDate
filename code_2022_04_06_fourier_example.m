% fourier series

myfun = @(x) x.^2 + 4*x -5; % y = x



x = -pi : 0.01 : pi;
N = 20;

val = an(myfun,0)/2;
for n = 1:N
    val = val + cos(n.*x) * an(myfun,n) + sin(n.*x) * bn(myfun,n);
end

close all;
figure; hold on;
plot(x, myfun(x));
plot(x, val);
legend('original function', 'approximated function');

% coefficient an
function coef = an(infun, n)
    coef = integral(@(x) infun(x).*cos(n*x), -pi, pi)/pi;
end

% coefficient bn
function coef = bn(infun, n)
    coef = integral(@(x) infun(x).*sin(n*x), -pi, pi)/pi;
end
