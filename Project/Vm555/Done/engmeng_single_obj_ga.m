clear
clc
tic;

J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction


ObjectiveFunction = @simple_fitness;
nvars = 5;    % Number of variables
LB = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001];   % Lower bound
UB = [1000 10 100 1 1]; % Upper bound

%Options
options = optimoptions(@ga);
options = optimoptions(options, 'PopulationSize', 50);
options = optimoptions(options, 'Generations', 25);
options = optimoptions(options, 'Display', 'iter');
%options = optimoptions(options, 'PlotFcns', {@gaplotbestindiv});
options = optimoptions(options,'MutationFcn',@mutationadaptfeasible,'PlotFcn',{@gaplotbestf}, ...
    'Display','iter');
    

% Next we run the GA solver.
[x,fval] = ga(@simple_fitness,nvars,[],[],[],[],LB,UB,[],options)



toc

function y = simple_fitness(x)
model = 'SBW_SystemModel';
% open_system(model);

in = Simulink.SimulationInput(model);

in = in.setVariable('x',x);
out = sim(in);
y = out.Error(end);
% Index2 = Consumption(end)
end