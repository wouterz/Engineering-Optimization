K_simple = [4500 1650 1100 2250 550 9300];
K = [K_simple(1)+K_simple(3)+K_simple(4), -K_simple(3), -K_simple(4);
    -K_simple(3), K_simple(2)+K_simple(3)+K_simple(5), -K_simple(5);
    -K_simple(4), -K_simple(5), K_simple(4)+K_simple(5)+K_simple(6)];
P = [1100 1800 3300];

opt = optimset();
[xopt, fopt] = fminsearch(@func, [0,0,0]', opt, K, P');

sprintf('PE %g - x1 %g, - x2 %g - x3 %g', fopt, xopt(1), xopt(2), xopt(3))

function y = func(x, K, P)
y = 0.5*x'*K*x - x'*P;
end