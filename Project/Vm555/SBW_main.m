clear
close all
clc

%% Model Parameters
J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

%% Design Parameters
% kp   = 200; % 0<=kp<=1000;
% ki   = 1;   % 0<=ki<=10;
% kd   = 20;  % 0<=kd<=100;
% kv   = 0;   % 0<=kv<=1;
% kdis = 0;   % 0<=kdis<=1;

kp   = 0.0001; % 0<=kp<=1000;
ki   = 0.0001;   % 0<=ki<=10;
kd   = 0.0001;  % 0<=kd<=100;
kv   = 0.0001;   % 0<=kv<=1;
kdis = 0.0001;   % 0<=kdis<=1;
x = [kp, ki, kd, kv, kdis];

%%
sim_model_name = 'SBW_SystemModel'


sim(sim_model_name);
Index1 = Error(end)
% Index2 = Consumption(end)


x = [200,1,20,0,0];

% open_system(sim_model_name)
% in = Simulink.SimulationInput(sim_model_name);
% in.setVariable('x', x2)
% in.applyToModel

% sim(in);

sim(sim_model_name);
% Index1 = Error(end)

Index1_2 = Error(end)
% Index2 = Consumption(end)

disp([Index1, Index1_2])

%%
figure()
plot(Angle(:,1),Angle(:,2)/pi*180,'r')
hold on
plot(Angle(:,1),Angle(:,3)/pi*180,'k')
xlabel('t(s)')
ylabel('angle(deg)')
legend('Reference','Tracking')
legend('boxoff')






