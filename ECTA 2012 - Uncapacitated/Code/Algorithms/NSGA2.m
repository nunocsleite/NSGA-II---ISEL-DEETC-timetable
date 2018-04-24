%-------------------------------------------------------------
%  Final Project of OMBAE PhD course, Winter semester 2011/12
%
%  Author: Nuno Miguel Leite, No. 59734/D
%-------------------------------------------------------------
function [NonDominated Rt] = NSGA2(Pt, Qt, Data, f, M, fLimits, crossProb, mutProb, ...
                                         reinsertionRate, cross, mutate, ...
                                         HCAlgorithm)
    % NSGA-II Algorithm
    %        
    generationNumb = 1;
    Rt = cell(Data.NIND*2,1);
    
    % Keep pareto fronts at iterations 1, 25, 50, 75, 100 and 125.
    NonDominated = cell(6,1);
    kk = 1;
    
    while (true)
        disp(sprintf('Generation: %d', generationNumb));

        % Step 1 - Combine parent and offspring populations 
        %     and create Rt = Pt U Qt. 
        for i = 1 : Data.NIND
            Rt{i} = Pt{i};
        end
        for i = 1 : Data.NIND
            Rt{Data.NIND+i} = Qt{i};
        end
              
        popSize = size(Rt,1);

        %     Perform a non-dominated sorting to Rt and identify
        %     different fronts: Fi, i = 1, 2, ..., etc.
        P = NonDominatedSorting(Rt, popSize, f{1}, f{2}); % TODO Passar array de handlers para função
        for i = 1 : length(P.Pk)-1
            P.Pk{i};
        end
        Fronts = P.Pk;
        % Each solution has the fields:
        %  1. A non-domination rank ri in the population.
        %  2. A local crowding distance (di) in the population.
        S.ri = zeros(1, popSize);
        S.di = zeros(1, popSize);
        
        % Step 2 - Set new population P(t+1) = empty. Set a counter i = 1.
        %     Until |P(t+1)| + |Fi| < N, perform P(t+1) = P(t+1) U Fi 
        %     and i = i+1.
        newPopIdxs = []; % P(t+1)
        i = 1;
        frontNumber = 1;
        while (i <= length(Fronts)-1 && length(newPopIdxs) + length(Fronts{i}) <= Data.NIND)
            newPopIdxs = [newPopIdxs Fronts{i}]; % P(t+1) = P(t+1) U Fi 
            % Actualize non-domination rank
            % Solution index
            m = Fronts{i};
            S.ri(m) = frontNumber;
            i = i+1;
            frontNumber = frontNumber+1;
        end
        % Last front Fi
        Fi = Fronts{i};
        % Actualize non-domination rank of the last front that can be
        % accommodated
        m = Fi;
        S.ri(m) = frontNumber;

        j = i+1;
        % Continue the ranking of the remainder fronts
        while (j <= length(Fronts)-1)
            % Actualize non-domination rank
            m = Fronts{j};
            S.ri(m) = frontNumber;
            j = j+1;
            frontNumber = frontNumber+1;
        end

        FrontIdxs = 1:i;
        % Step 3 - Perform the Crowding-sort(Fi, <c) procedure and include
        %     the most widely spread (N - |P(t+1)|) solutions by using the
        %     crowding distance values in the sorted Fi to P(t+1).
        %
        % Computes Crowding distances
        S = CrowdingDistances(Rt, S, Fronts, FrontIdxs, f, M, fLimits);
        % Crowding sort
        CrowDist = S.di;
        FiSorted = CrowdingSort(Fi, CrowDist);
        % Number of remainder solutions of Fi to include in the population
        K = Data.NIND-length(newPopIdxs);
        % Include the most widely spread (N - |P(t+1)|) solutions by using 
        % the crowding distance values in the sorted Fi to P(t+1).
        for k = 1 : K
           newPopIdxs = [newPopIdxs FiSorted(k)]; 
        end
 
        if (generationNumb >= 2)
            % Improve solutions using local search
            Rt = improvePopulation(S, Data, Rt, HCAlgorithm);
        end   
             
        if (generationNumb == 1)
            NonDominated{kk} = nonDominated(Rt, newPopIdxs, S);
            kk = kk+1;
        elseif (mod(generationNumb, 25) == 0)
            NonDominated{kk} = nonDominated(Rt, newPopIdxs, S);
            kk = kk+1;
        end
    
        Pt = Rt(newPopIdxs, :);
      
        if (generationNumb == Data.MAXGEN)
            break;
        end
        % Step 4 - Create offspring population Q(t+1) from P(t+1) by 
        %      using the crowded tournament selection, crossover
        %      and mutation operators.
        %
        % Crowded Tournament Selection
        offspringIdxs = select(newPopIdxs, S);
        % build Qt
        Qt = cell(Data.NIND, 1);
        for i = 1 : Data.NIND
            Qt{i} = Rt{offspringIdxs(i)};
        end
        % Cross P(t+1)
        Qt = cross(Data, Qt, crossProb);
        % Mutate some solution from P(t)
        Qt = mutate(Data, Qt, mutProb, reinsertionRate);

        generationNumb = generationNumb+1;
    end
end

function pop = improvePopulation(S, Data, pop, HCAlgorithm)
    % The local search operator is applied to chromosomes
    % selected based on a tournament selection scheme,
    % where all the chromosomes in the population are randomly
    % grouped into fours and from each group, the chromosome
    % with the smallest rank is selected. Only a quarter of the
    % population will undergo local exploitation.    
    popSize = length(pop);
    p = randperm(popSize);
    
    for i = 1 : 4 : popSize
        ranks = [S.ri(i), S.ri(i+1), S.ri(i+2), S.ri(i+3)];
        [M, I] = min(ranks);
        improvIdx = i+I-1;
        % Improve with Local search
        pop{improvIdx} = HCAlgorithm(pop{improvIdx}, Data);
    end
end
%////////////////////////////////////////////////////////////


% Get non-dominated solutions
function Pt = nonDominated(Rt, newPopIdxs, S)
    nonDominatedIdxs = find(S.ri(newPopIdxs) == 1); % Find rank 1 solutions
    for i = 1 : length(nonDominatedIdxs)
        Pt{i} = Rt{newPopIdxs(nonDominatedIdxs(i))};
    end
end

% Crowding-sort(Fi, <c) - Crowding Distance Assignment Procedure
function S = CrowdingDistances(Rt, S, Fronts, FrontIdxs, f, M, fLimits)
    % Compute crowding distances for solutions in each front
    for i = 1 : length(FrontIdxs)
        % Step C1 - Call the number of solutions in F as l = |F|.
        %      For each i in the set, first assign di=0.
        Fi = Fronts{i};
        l = length(Fi);
        % Step C2 - For each objective function m = 1,2,...,M, sort the
        %      set in worse order of fm or, find the sorted indices
        %      vector: I^m = sort(fm, >).
        I = cell(1,M);

        for m = 1 : M
            fmValues = zeros(l,1);
            for z = 1 : l
                fmValues(z) = f{m}(Rt{Fi(z)});
            end
            [sortedVals, Idxs] = sort(fmValues);
            I{m} = Fi(Idxs);
            % Crowding distance calculation
            firstIdx = I{m}(1);
            lastIdx = I{m}(l);
            fmMin = fLimits(2*m-1);
            fmMax = fLimits(2*m);
            den = fmMax-fmMin;
            S.di([firstIdx lastIdx]) = inf;
            for k = 2 : l-1
                s = I{m}(k);
                S.di(s) = S.di(s) + (sortedVals(k+1)-sortedVals(k-1))/den;
            end
        end
     end
end

% Crowding-sort(Fi, <c) procedure.
% Sorts Fi acording to crowding distance values.
function FiSorted = CrowdingSort(Fi, CrowDist)
    [sortedVals, Idxs] = sort(CrowDist(Fi), 'descend');
    FiSorted = Fi(Idxs);
end

function offspringIdxs = select(newPopIdxs, S)
    popSize = length(newPopIdxs);
    popRandIdxs = randi([1, popSize], 1, popSize);
    % Binary Crowded tournament selection
    for i = 1 : popSize
        Si = newPopIdxs(i);
        Sj = newPopIdxs(popRandIdxs(i));
        Sol = CrowdedTournament(S, Si, Sj);
        offspringIdxs(i) = Sol;
    end 
end

%
% Crowded Tournament Selection Operator
%
% The crowded comparison operator (<c) compares two solutions and
% returns the winner of the tournament. It assumes that every 
% solution i has two attributes:
%    1. A non-domination rank ri in the population.
%    2. A local crowding distance (di) in the population.
%
% The crowding distance di of a solution i is a measure of the search
% space around i which is not occupied by any other solution in
% the population. 
%
% Crowded Tournament Selection Operator: A solution i wins a tournament
% with another solution j if any of the following conditions are true:
%    1. If solution i has a better rank, that is, ri < rj.
%    2. If they have the same rank but solution i has a better crowding
%       distance than solution j, that is, ri == rj and di > dj.
function Sol = CrowdedTournament(S, i, j)
    if (S.ri(i) < S.ri(j)|| ...
        (S.ri(i) == S.ri(j) && S.di(i) > S.di(j)))
        Sol = i;
    else
        Sol = j;
    end
end

