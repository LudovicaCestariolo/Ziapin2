function [VOI,STATES,ALGEBRAIC,CONSTANTS] = main()
    
% clear; close all; clc;

protocol = 2;   % protocol = 0 --> Control
                % protocol = 1 --> ZiaPina2 electric stimulation
                % protocol = 2 --> ZiaPina2 dark-light 20 ms
                % protocol = 3 --> ZiaPina2 dark-light 200 ms

[VOI,STATES,ALGEBRAIC,CONSTANTS] = solveModel(protocol);
end

function [VOI,STATES,ALGEBRAIC,CONSTANTS] = solveModel(protocol)

% Create ALGEBRAIC of correct size
global algebraicVariableCount;  
algebraicVariableCount = 109;

% Initialise constants and state variables
[INIT_STATES,CONSTANTS] = initConsts(protocol);

% Set timespan to solve over
tspan = [0,1000];

% Set numerical accuracy options for ODE solver
options = odeset('RelTol',1e-02,'AbsTol',1e-02,'MaxStep',1);

% Solve model with ODE solver
[VOI,STATES] = ode15s(@(VOI,STATES)computeRates(VOI,STATES,CONSTANTS,protocol),tspan,INIT_STATES,options);

% Compute algebraic variables
[RATES,ALGEBRAIC] = computeRates(VOI,STATES,CONSTANTS,protocol);
ALGEBRAIC = computeAlgebraic(ALGEBRAIC,CONSTANTS,STATES,VOI,protocol);

% Plot state variables against variable of integration
figure(1);
plot(VOI, STATES(:,1));

end