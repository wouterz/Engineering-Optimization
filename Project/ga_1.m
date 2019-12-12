clear
% clc

ObjectiveFunction = @simple_fitness;
nvars = 2;    % Number of variables
LB = [0.00001 0.00001];   % Lower bound
UB = [inf inf];  % Upper bound

%Options
% options = optimoptions;
options = optimoptions('ga', 'PopulationSize', 50);
options = optimoptions(options, 'Generations', 200);
options = optimoptions(options, 'Display', 'final');
options = optimoptions(options, 'PlotFcns', {@gaplotbestindiv,
@gaplotscorediversity});
options = optimoptions(options,'MutationFcn',@mutationadaptfeasible,'PlotFcn',{@gaplotbestf,@gaplotmaxconstr}, ...
    'Display','iter');
    
% X0 = [0 0]; % Start point (row vector)
% options.InitialPopulationMatrix = X0;
% Next we run the GA solver.
[x,fval] = ga(@simple_fitness,nvars,[],[],[],[],LB,UB, ...
    @simple_constraint,options)

function y = simple_fitness(x)
% m1 = x(1);
% m2 = x(2);
% [T,X,Y] = sim('test1.mdl');

    y = (2*x(1)+3*x(2)-x(1)^3-2*x(2)^2);
end
function [c, ceq] = simple_constraint(x)
c(1) = x(1)+3*x(2)-6;
c(2) = 5*x(1)+2*x(2)-10;
  
   %Equality constraints
   ceq = [];
end
