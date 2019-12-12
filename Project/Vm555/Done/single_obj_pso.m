tic;
%Single Objective PSO
%Phase 1 - Remove all previous values from the workspace
clc;
clear;
close all;

%Declare constants for the simulink model here
J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

%Phase 2 - Define the problem
CostFunction = @single_Obj;
%Define the number of decision variables
nVar = 5;
%Define decision matrix size
VarSize=[1 nVar];
%Define the normalized lower and upper bound of the problem
VarMin=0;
VarMax=1;

%Phase 3 - Define the PSO parameters
%Define max number of iterations
MaxIt=25; 
%Define population size
nPop = 50;
%Intertia coefficient
w = 0.5;
%Damping coefficient
wdamp = 0.99;
%Personal learning coefficient
c1 = 1;
%Global learning coefficient
c2 = 1;

%Phase 4 - Initialize the particles and the information
empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

%Create population array for all the particle information
particle=repmat(empty_particle,nPop,1);

%Initialize the global best
GlobalBestCost = inf;
GlobalBest.Cost = inf;

%Initialize population members
for i=1:nPop

    %Generate Random solution
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    %Initialize velocity
    particle(i).Velocity=zeros(VarSize);
    %Apply small velocity to kickstart 
    particle(i).Velocity = 0.1;
    
    %Evaulate cost value of solution
    particle(i).Cost=CostFunction(particle(i).Position);
    
    %Update personal best and cost to current location
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    %Update global best
    if particle(i).Best.Cost < GlobalBest.Cost
        GlobalBest = particle(i).Best;
    end 

end

%Get best cost value for each iteration
BestCosts = zeros(MaxIt,1);

%Phase 5 - Main loop of PSO

for it = 1:MaxIt
    
    for i = 1:nPop
       %Determine the updated velocity of the particles
        particle(i).Velocity = w*particle(i).Velocity ...
                + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
                + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);
       
       %Determine the updated position of the particles
%        particle(i)
       particle(i).Position = particle(i).Position + particle(i).Velocity;
       
       % Apply Lower and Upper Bound Limits
       particle(i).Position = min(particle(i).Position, VarMax);
       particle(i).Position = max(particle(i).Position, VarMin);
       
       %Determine the updated cost function
       particle(i).Cost = CostFunction(particle(i).Position);

       %Update the personal best
       if particle(i).Cost < particle(i).Best.Cost
           particle(i).Best.Position = particle(i).Position;
           particle(i).Best.Cost = particle(i).Cost;
           
           %Update global best
           if particle(i).Best.Cost < GlobalBest.Cost
               GlobalBest = particle(i).Best;
           end 
       end 
       
       
    end 
    
    %Store best cost value 
    BestCosts(it) = GlobalBest.Cost;
%     disp(['asdasd', num2str(GlobalBest.Cost), BestCosts(it)])

    
    %Display iteration info
    disp(['Iteration' num2str(it) ' Best Cost =' num2str(BestCosts(it))]);
    fprintf('Params: %g %g %g %g %g\n', GlobalBest.Position.*[1000, 10, 100,1,1])
    w = w * wdamp;
end 

    
%Phase 7 - Display graphs and results
figure
plot(BestCosts);
xlabel('Iterations');
ylabel('Cost Function');

toc
%%Define the function
function z = single_Obj(x)
model = 'SBW_SystemModel';
J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction
x = abs(x) .*[1000, 10, 100,1,1];

in = Simulink.SimulationInput(model);
in = in.setVariable('x',x);
out = sim(in);
z = out.Error(end);
%f1 = z;
    %z = [f1];
end

