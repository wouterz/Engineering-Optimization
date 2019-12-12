clear
close all
clc

%% Model Parameters
J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

%% Design Parameters
kp   = 200; % 0<=kp<=1000;
ki   = 1;   % 0<=ki<=10;
kd   = 20;  % 0<=kd<=100;
kv   = 0;   % 0<=kv<=1;
kdis = 0;   % 0<=kdis<=1;

%%
sim('SBW_SystemModel');
Index1 = Error(end)+ max(trackerror);
Index2 = Consumption(end);

%%
figure()
plot(Angle(:,1),Angle(:,2)/pi*180,'r')
hold on
plot(Angle(:,1),Angle(:,3)/pi*180,'k')
xlabel('t(s)')
ylabel('angle(deg)')
legend('Reference','Tracking')
legend('boxoff')





