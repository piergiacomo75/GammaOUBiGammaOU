classdef GOUBGOU
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = public)
        lambdaPoisson
        lambdaJumps
        kappaOU
        x_0
    end
    
    methods
        function obj = GOUBGOU(lambdaPoisson, lambdaJumps, kappaOU, x_0)

            if nargin < 4
                x_0 = 0;
            end
            obj.lambdaPoisson = lambdaPoisson;
            obj.lambdaJumps = lambdaJumps;
            obj.kappaOU = kappaOU;
            obj.x_0 = x_0;
            
        end
        
        function trajectory = SimulateTrajectory(obj, timeGrid, numberOfSimulations, simulationMethodStr, numberOfCoefficients)
            
            if nargin < 5
                numberOfCoefficients = 0;
            end
            
            deltaGrid = diff(timeGrid);
            numberOfTimeSteps = numel(deltaGrid);
            
            simulationMethod = ResolveIncrementMethod(simulationMethodStr, numberOfCoefficients);
            trajectory = zeros(numberOfTimeSteps, numberOfSimulations);
            trajectory(1, :) = obj.x_0;
            for t = 1:numberOfTimeSteps
                increment = simulationMethod.GenerateIncrement(obj.lambdaPoisson, obj.lambdaJumps, obj.kappaOU, deltaGrid(t), numberOfSimulations);                
                trajectory(t+1, :) = trajectory(t, :) * exp(-obj.kappaOU * deltaGrid(t)) + increment'; 
            end
            
        end
    end
end

