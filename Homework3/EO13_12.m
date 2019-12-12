A = [1,1,1];
B = [-1, 2, 4];
C = [2, 3, 4];
D = [-3,-4,1];
points = [A; B; C; D];

X0 = [4,4,4,4];
%x y z r
lb = [-inf, -inf, -inf, 0];
ub = [];
options = optimset('Algorithm', 'active-set');

[x,f] = fmincon(@NonLinF, X0, [], [], [], [], ...
    lb, ub, @NonLinConst, options, points, size(points));

sprintf('%g - %g - %g  - %g', x)
function y = NonLinF(x, pts, rad)
y = 2*x(4);
end

function [C, Ceq] = NonLinConst(x, pts, rad)
for i = 1:rad
    C(i) = (x(1)-pts(i,1))^2 ...
            +(x(2)-pts(i,2))^2 ...
            +(x(3)-pts(i,3))^2 ...
            -x(4)^2;
end
Ceq = [];
end
