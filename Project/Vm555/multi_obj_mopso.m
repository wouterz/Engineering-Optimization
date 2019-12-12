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
MaxIt=250; 
%Define population size
nPop = 5;
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
empty_particle.Best.Cost2=[];
empty_particle.Delete=0;

%Create population array for all the particle information
particle=repmat(empty_particle,nPop,1);

%Initialize the global best
GlobalBest.Cost = inf;
GlobalBest.Cost2 = inf;

rep = [];

%Initialize population members
for i=1:nPop

    %Generate Random solution
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    %Initialize velocity
    particle(i).Velocity=zeros(VarSize);
    %Apply small velocity to kickstart 
    particle(i).Velocity = 0.1;
    
    %Evaulate cost value of solution
    c = CostFunction(particle(i).Position);
    particle(i).Cost = c(1);
    particle(i).Cost2 = c(2);
    
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
       particle(i).Position = particle(i).Position + particle(i).Velocity;
       
       % Apply Lower and Upper Bound Limits
       particle(i).Position = min(particle(i).Position, VarMax);
       particle(i).Position = max(particle(i).Position, VarMin);
       
       %Determine the updated cost function
       c = CostFunction(particle(i).Position);
       particle(i).Cost = c(1);
       particle(i).Cost2 = c(2);
               
       rep = [rep particle(i)];

       if numel(rep) > 1
       
            rep = recheckDominance(rep);
            fprintf('rep size %i -----  ', numel(rep))
            rep=rep(~[rep.Delete]);
            fprintf('rep size after delete  %i \n', numel(rep))

       end
       

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
    
    w = w * wdamp;
end 

    
%Phase 7 - Display graphs and results
figure
plot(BestCosts);
xlabel('Iterations');
ylabel('Cost Function');


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
z(1) = out.Error(end);
z(2) = out.Consumption(end);

%f1 = z;
    %z = [f1];
end

function y = dominates(point1, point2)
    y = false;

    if point1.Cost < point2.Cost ...
       && point1.Cost2 < point2.Cost2
        y = true;
    end
end

function y = recheckDominance(current_leaders)
    n=numel(current_leaders);
    for i=1:n-1
        for j=i+1:n
            
            if dominates(current_leaders(i),current_leaders(j))
               current_leaders(j).Delete = 1;
            end
            
            if dominates(current_leaders(j),current_leaders(i))
               current_leaders(i).Delete = 1;
            end
            
        end
    end
    y = current_leaders;
end

