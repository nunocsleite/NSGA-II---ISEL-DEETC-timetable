%-------------------------------------------------------------
%  Final Project of OMBAE PhD course, Winter semester 2011/12
%
%  Author: Nuno Miguel Leite, No. 59734/D
%-------------------------------------------------------------

% An O(M N^2) Non-Dominated Sorting Algorithm
function P = NonDominatedSorting(Pop, popSize, f1, f2) 
    % ALGORITHM:
    %
    % First, for each solution we calculate two entities: (1) domination
    % count ni, the number of solutions which dominate solution i, 
    % and (2) Si, a set of solutions which the solution i dominates.
    % This requires O (M N^2) comparisons.
    % At the end of this procedure, all solutions in the first 
    % non-dominated front will have their domination count as zero.
    %   
    % Now, for each of these solutions (each solution i with ni = 0),    
    % we visit each member (j) of its set Si and reduce its domination 
    % count by one. In doing so, if for any member j the domination count
    % becomes zero, we put it in a separate list P'. After such
    % modifications on Si are performed for each i with ni = 0, all
    % solutions of P' would belong to the second non-dominated front.
    % The above procedure can be continued with each member of P' and the
    % third non-dominated front can be identified. This process continues
    % until all solutions are classified.
    %
    % STEPS OF THE ALGORITHM:
    %
    % Step 1 - For each i in P, ni = 0 and initialize Si = empty.
    % For all j <> i and i in P, perform Step 2 and then 
    % proceed to Step 3.

    % Step 2 - If i <= j, update Si = Si U {j}. Otherwise, 
    % if j <= i, set ni = ni + 1.
    
    % Step 3 - If ni = 0, keep i in the first non-dominated front P1. 
    % Set front counter k = 1.
    
    % Step 4 - While Pk <> 0, perform the following steps.
    
    % Step 5 - Initialize Q = empty for storing next non-dominated
    % solutions. For each i in Pk and for each j in Si,
    %   Step 5a - Update nj = nj - 1 
    %   Step 5b - If nj == 0, add j in Q, or perform Q = Q U {j}
    
    % Step 6 - Set k = k + 1 and Pk = Q. Go to Step 4. 
    % -------------------------------------------------------------------
    
    % Step 1 - For each i in P, ni = 0 and initialize Si = empty.
    P.ni = zeros(1, popSize); % ni, is the domination count
    P.Si = cell(1, popSize);  % Si, set of solutions which the 
                              % solution i dominates
    % For all j <> i and i in P, perform Step 2 and then 
    % proceed to Step 3.
    for i = 1 : popSize
        for j = 1 : popSize
            if (j ~= i)
                % Step 2 - If i <= j, update Si = Si U {j}. Otherwise, 
                % if j <= i, set ni = ni + 1.
                if (dominates(Pop, i, j, f1, f2))    
                    P.Si{i} = [P.Si{i} j];
                else if (dominates(Pop, j, i, f1, f2))    
                    P.ni(i) = P.ni(i) + 1;
                    end
                end
            end
        end
    end
        
    % Step 3 - For all i, If ni = 0, keep i in the 
    %          first non-dominated front P1. 
    idx = find(P.ni == 0);
    P.Pk{1} = idx; % Front 1
    % Set front counter k = 1.
    k = 1;    
    % Step 4 - While Pk not empty, perform the following steps.
    while (~isempty(P.Pk{k}))
        % Step 5 - Initialize Q = empty for storing next non-dominated
        % solutions. 
        Q = [];
        % For each i in Pk and for each j in Si,
        for ii = 1 : length(P.Pk{k})
            i = P.Pk{k}(ii);
            for jj = 1 : length(P.Si{i})
                j = P.Si{i}(jj);
                %   Step 5a - Update nj = nj - 1 
                P.ni(j) = P.ni(j) - 1;
                %   Step 5b - If nj == 0, add j in Q, or perform Q = Q U {j}
                if (P.ni(j) == 0)
                    Q = [Q j];
                end
            end
        end
        % Step 6 - Set k = k + 1 and Pk = Q. Go to Step 4. 
        k = k + 1;
        P.Pk{k} = Q;
    end
end

function b = dominates(Pop, i, j, f1, f2) % TODO: tornar genérico
    xi = Pop{i};
    xj = Pop{j};
    if (notWorse(xi, xj, f1, f2) && better(xi, xj, f1, f2) )
        b = true;
    else
        b = false;
    end
end

function b = notWorse(xi, xj, f1, f2) % TODO: tornar genérico
    if (f1(xi) <= f1(xj) && f2(xi) <= f2(xj))
        b = true;
    else 
        b = false;
    end
end

function b = better(xi, xj, f1, f2) % TODO: tornar genérico
    if (f1(xi) < f1(xj) || f2(xi) < f2(xj))
        b = true;
    else 
        b = false;
    end
end
