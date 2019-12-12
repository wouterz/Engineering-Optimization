clc;
clear;
close all;

tic;

J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

%% Problem Definition

CostFunction= @multi_obj_func;      % Cost Function

nVar=5;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax=1;          % Upper Bound of Variables


%% MOPSO Parameters

MaxIt=25;           % Maximum Number of Iterations

nPop=50;            % Population Size

nRep=100;            % Repository Size

w=0.5;              % Inertia Weight
wdamp=0.99;         % Intertia Weight Damping Rate
c1=1;               % Personal Learning Coefficient
c2=2;               % Global Learning Coefficient

alpha=0.1;          % Inflation Rate

beta=2;             % Leader Selection Pressure
gamma=2;            % Deletion Selection Pressure

mu=0.1;             % Mutation Rate

%% Initialization

empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.IsDominated=false;

particle=repmat(empty_particle,nPop,1);
i = 1;
while i <= nPop
    % Generate a Random Solution
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    % Initalize Velocity
    particle(i).Velocity=zeros(VarSize);
    % Evaulation
    particle(i).Cost=CostFunction(particle(i).Position);
    
    % Update Best to its current initial location
    particle(i).Best.Position=particle(i).Position;
    % Update best cost
    particle(i).Best.Cost=particle(i).Cost;
    
    % Forice initial population to have valid cost
    if particle(i).Cost(1) < 9
        i= i+1;
    end
    
end

% Evaluate domination
particle=evaluateDomination(particle);

best_particles=particle(~[particle.IsDominated]);

%% MOPSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        leader=best_particles(randi([1 numel(best_particles)]));
        
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-particle(i).Position);
        
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        particle(i).Position = max(particle(i).Position, VarMin);
        particle(i).Position = min(particle(i).Position, VarMax);
        
        particle(i).Cost = CostFunction(particle(i).Position);
        
        if dominates(particle(i),particle(i).Best)
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
        elseif dominates(particle(i).Best,particle(i))
            % Nothing
            
        else
            if rand<0.5
                particle(i).Best.Position=particle(i).Position;
                particle(i).Best.Cost=particle(i).Cost;
            end
        end
        
    end
    
    % Add particles that are not dominated to the repository
    best_particles=[best_particles
         particle(~[particle.IsDominated])];
    
    % Evaluate if any particles are dominated
    best_particles=evaluateDomination(best_particles);
    
    % Drop dominated particles after evaluation
    best_particles=best_particles(~[best_particles.IsDominated]);
    
    % Remove random particle if best_particles is full
    while numel(best_particles)>nRep
        best_particles(randi([1 numel(best_particles)]))=[];        
    end
    
    % Plot Costs
    figure(1);
    xlabel('Error');
    ylabel('Consumption');
    %Plot pareto front points
    rep_costs=[best_particles.Cost];
    plot(rep_costs(1,:),rep_costs(2,:),'r*');
    grid on
    
    % Display information
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(best_particles))]);
    
    % dampening
    w=w*wdamp;
    
end

toc

function y = multi_obj_func(x)
model = 'SBW_SystemModel';

J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

in = Simulink.SimulationInput(model);

in = in.setVariable('x',x);

out = sim(in);
f1 = out.Error(end);
f2 = out.Consumption(end);
y = [f1
    f2];
end

function b = dominates(x,y)
    if isstruct(x)
        x=x.Cost;
    end
    if isstruct(y)
        y=y.Cost;
    end
    % Atleast one of the inequalities must be strict, others can be equal
    % or smaller
    b = any(x<y) && all(x<=y);

end

function population = evaluateDomination(population)
    n=numel(population);
    % Reset all dominated values to false before beginning
    for i=1:n
        population(i).IsDominated=false;
    end
    
    % Pairwise evaluation of each combination
    for i=1:n-1
        for j=i+1:n
            if dominates(population(i),population(j))
               population(j).IsDominated=true;
            end
            
            if dominates(population(j),population(i))
               population(i).IsDominated=true;
            end
            
        end
    end
end