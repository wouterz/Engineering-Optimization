function y = singleObjFunc(x)
model = 'SBW_SystemModel';
% model = 'test2';
J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

in = Simulink.SimulationInput(model);
in = in.setVariable('x',x);

out = sim(in);
% y = out.simout(end);
y = out.Error(end);
% Index2 = Consumption(end)
end

        