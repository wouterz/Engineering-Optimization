E = 1724;
L = 127;
t = 25.4;

R = [27.9 43.2 53.3 73.7 99.1];
F = [385 536 2317 7820 18230];


x = lsqcurvefit(@func, [1 1], R, F, [], [], [], E, L, t);
sprintf('a %g - b %g', x(1), x(2))

figure;
hold on;

plot(R,F, 'o')
plot(R, func(x, R,E,L,t))
function F = func(x, R, E, L, t)
F = (pi^x(1)*E*R.^x(2)*t^(4-x(2)))/(4*L^2);
end