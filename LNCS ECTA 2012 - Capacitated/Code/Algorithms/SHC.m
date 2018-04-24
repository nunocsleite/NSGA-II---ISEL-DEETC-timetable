%-------------------------------------------------------------
%  Final Project of OMBAE PhD course, Winter semester 2011/12
%
%  Author: Nuno Miguel Leite, No. 59734/D
%-------------------------------------------------------------
function Res = SHC(tmax, T, data, initialSolution, getRandomNeigh, ...
    evalFunc)
    % SHC - Stochastic hill climber (minimization of f) 
    %  
    % Choose a solution u (at random) and compute fu = f(u) 
    % Repeat step 1 to 3 (until t = tmax) 
    %   Step 1  Find (at random) a neighbor of u, say v. 
    %   Step 2  Evaluate v, i.e., compute fv = f(v) 
    %   Step 3  Make u = v with probability p(fu,fv), see Eq.(1); 
    %   Make t = t+1. 
    % 
    sense = 'minimize';
    
	% Number of evaluations
    numEvaluations = 0;
    % Choose a solution u (at random) and compute fu = f(u) 
    u = initialSolution;
    fu = evalFunc(u, data);
    % Increment number of evaluations
    numEvaluations = numEvaluations + 1;
    z = 1;
    F(z) = fu;
    z = z+1; 
    orig = fu;
    
    % Repeat step 1 to 3 (until t = tmax) 
    t = 1;
    while (t <= tmax)
        % Step 1  Find (at random) a neighbor of u, say v. 
        v = getRandomNeigh(data, u);
        % Step 2 Evaluate v, i.e., compute fv = f(v) 
	    fv = evalFunc(v, data);
	    % Increment number of evaluations
        numEvaluations = numEvaluations + 1;
        
        % Step 4 Make u = v with probability p(fu,fv)
        prob = p(fu, fv, T, sense);
        
        x = myRand();
        if (x < prob)
            % Accept this solution
            u = v;
            fu = fv;
        end
        F(z) = fu;
        z = z+1;
        % Make t = t+1. 
        t = t + 1;
    end

%     orig
%     fu

    Res = struct('T', T, 'tmax', tmax, 'NumEvaluations', numEvaluations, ...
        'Cost', fu, 's', u, 'u', u);
end


%/////////////////////////////////////////////////////

% Random number generator
function r = myRand()
    r = rand;
end

% Probability function
function x = p(fu, fv, T, sense)
    % Minimization problem
    %x = 1 / (1 + exp((fv-fu) / (T*fu)));
    % Maximization problem
    %x = 1 / (1 + exp((fu-fv) / (T*fu)));
    if strcmp(sense, 'maximize')
        factor = -1;
    elseif strcmp(sense, 'minimize')
        factor = 1;
    end
    x = 1 / (1 + exp(factor*(fv-fu) / (T*fu)));
end

