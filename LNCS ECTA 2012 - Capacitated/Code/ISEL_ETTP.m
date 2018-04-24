%-------------------------------------------------------------
%  Final Project of OMBAE PhD course, Winter semester 2011/12
%
%  Author: Nuno Miguel Leite, No. 59734/D
%-------------------------------------------------------------
function ISEL_ETTP
    % ISEL Examination Timetabling Problem
    % 
    clc; % Clear screen
    %
    % ISEL/DEETC department test data
    %
    Data = importfile('Relacoes UC 0910i.xls');
    
    %Data.Classes
    %size(Data.Classes)
    %pause
    
    
    ExamsToRemove = {'LIC', 'SEAD2', 'AC2', 'FIA', 'PCM', 'CIA', 'PSTR', ...
        'SEAS', 'CEE', 'PRC', 'SEADI', 'RM', 'STDS'}; % PRC does not have exam
    Data = removeExamsFromData(Data, ExamsToRemove);
    
    ExamCounts = computeNumberStudents(Data.ConflictMatrix);
    Data.ExamCounts = ExamCounts;
   
    tic; % Time measurement

    % Periods range definition. There are 3 period in each day
    % and one period at saturday.
    NumberWeeks = 3; % 1st epoch
    periods = NumberWeeks*6; % 1 period per day
    % Min number of periods.
%     Data.MinPeriods = periods;
    Data.MinPeriods = periods-4;

    % Max number of periods.
%     Data.MaxPeriods = periods;
    Data.MaxPeriods = periods+4;
    % Period range
    Data.Range = [Data.MinPeriods Data.MaxPeriods];
    % Population size
    Data.NIND = 40;
    %Data.NIND = 2;
    % Number of generations
    Data.MAXGEN = 125;
    % Total Period seating capacity    
    Data.S = 600;
    % Objective functions
    f = cell(1,2);
    f{1} = @f1;
    f{2} = @f2;
    % Number of objectives
    M = length(f);
    % Objective function's limits
    fLimits = [0 inf 0 inf];
    
    crossProb = 1;
    mutProb = 0.2;
%     mutProb = 1;
    % Rate at which exams are reinserted for each chromosome 
    % selected on mutation rate.
    reinsertionRate = 0.02;
    
    % Create rooms
    Data.Rooms = createRooms();
    
    %printRooms(Data.Rooms);
    
    
    % Generate initial populations Pt and Qt
    Pt = initialPopulation(Data);
    Qt = initialPopulation(Data);
     
    [NonDominated Rt] = NSGA2(Pt, Qt, Data, f, M, fLimits, crossProb, mutProb, ...
        reinsertionRate, @cross, @mutate, @HCAlgorithm);
     
    toc % Finish time measurement
    
    printPopulation(Data, Pt);
    printPopulation(Data, Qt);
    
    FinalFront = NonDominated{6};
    for i = 1 : length(FinalFront)
        t = computeNumClashes(Data, FinalFront{i});
        if (FinalFront{i}.NumClashes ~= t.NumClashes)
            disp('Error')
            t.NumClashes
            FinalFront{i}.NumClashes
            t.NumPeriods
            FinalFront{i}.NumPeriods
            pause
        end
    end
    
    showPareto(NonDominated, Rt);
    
    printFinalPopulation(Data, NonDominated{6});
    printExcelTimetables(Data, NonDominated{6});

    FinalPop = NonDominated{6};
    for i = 1 : length(FinalPop)
        % Assign rooms
        timetable = FinalPop{i};
        if (timetable.NumPeriods == 18)
            fprintf('# periods = %d', timetable.NumPeriods);
            assignRooms(Data, timetable);
        end
    end
    
end

function showPareto(NonDominated, Rt)
    close all;
    for i = 1 : 6
        % Current pareto front
        Fi = NonDominated{i};
        x = zeros(length(Fi),1);
        y = zeros(length(Fi),1);
        for j = 1 : length(Fi)
            s = Fi{j};
            x(j) = s.NumPeriods;
            y(j) = s.NumClashes;
        end
        subplot(3,2,i),
        plot(x, y, '+');
        for j = 1 : length(Fi)
            s = sprintf('(%d,%d)', x(j), y(j)); 
            text(x(j)-.5, y(j)+80, s);
        end
        axis([12, 24, 0, 1500]);
    end
end


function printExcelTimetables(Data, Pop)
    disp('Generating excel timetables')
    for i = 1 : length(Pop)
        disp(sprintf('Solution #%d', i));
        x = Pop{i};
%         fprintf('# periods = %d\n', x.NumPeriods);
%         fprintf('# Clashes = %d\n', x.NumClashes);
        filename = strcat('DEETC_', num2str(x.NumPeriods), '_', num2str(x.NumClashes));
        generateExcelTimetables(filename, Data, x);
    end
end


%////////////////////////////////////////////////////////////////////////
%                           Rooms allocation
%////////////////////////////////////////////////////////////////////////


function Rooms = createRooms()
    RoomSize = 30;
    Rooms = cell(1, RoomSize);
    % A.2.03	&	50
    Rooms{1}.designation = 'A.2.03';
    Rooms{1}.capacity =	50;
    Rooms{1}.Free = true;
    % A.2.08	&	40
    Rooms{2}.designation = 'A.2.08-A.2.09';
    Rooms{2}.capacity =	40+40;
    Rooms{2}.Free = true;
    % A.2.09	&	40
%    Rooms{1}.designation = 'A.2.09';
%    Rooms{1}.capacity =	40;
    % A.2.10	&	40
    Rooms{3}.designation = 'A.2.10-A.2.11';
    Rooms{3}.capacity =	40+40;
    Rooms{3}.Free = true;
    % A.2.11	&	40
%    Rooms{1}.designation = 'A.2.11';
%    Rooms{1}.capacity =	40;
    % A.2.12	&	40
    Rooms{4}.designation = 'A.2.12-A.2.13';
    Rooms{4}.capacity =	40+40;
    Rooms{4}.Free = true;
    % A.2.13	&	40	
%    Rooms{1}.designation = 'A.2.13';
%    Rooms{1}.capacity =	40;
    % A.2.16	&	45
    Rooms{5}.designation = 'A.2.16-A.2.18';
    Rooms{5}.capacity =	45+45;
    Rooms{5}.Free = true;
    % A.2.18	&	45
%    Rooms{1}.designation = 'A.2.18';
%    Rooms{1}.capacity =	45;
    % C.2.14	&	47
    Rooms{6}.designation = 'C.2.14';
    Rooms{6}.capacity =	47;
    Rooms{6}.Free = true;
    % C.2.21	&	16
    Rooms{7}.designation = 'C.2.21';
    Rooms{7}.capacity =	16;
    Rooms{7}.Free = true;
    % C.2.22	&	47	
    Rooms{8}.designation = 'C.2.22';
    Rooms{8}.capacity =	47;
    Rooms{8}.Free = true;
    % C.2.23	&	48
    Rooms{9}.designation = 'C.2.23';
    Rooms{9}.capacity =	48;
    Rooms{9}.Free = true;
    % C.3.07	&	75
    Rooms{10}.designation = 'C.3.07';
    Rooms{10}.capacity = 75;
    Rooms{10}.Free = true;
    % C.3.14	&	36	
    Rooms{11}.designation = 'C.3.14';
    Rooms{11}.capacity = 36;
    Rooms{11}.Free = true;
    % C.3.15	&	40	
    Rooms{12}.designation = 'C.3.15';
    Rooms{12}.capacity = 40;
    Rooms{12}.Free = true;
    % C.3.16	&	40	
    Rooms{13}.designation = 'C.3.16';
    Rooms{13}.capacity = 40;
    Rooms{13}.Free = true;
    % G.0.08	&	30
    Rooms{14}.designation = 'G.0.08';
    Rooms{14}.capacity = 30;
    Rooms{14}.Free = true;
    % G.0.13	&	30
    Rooms{15}.designation = 'G.0.13';
    Rooms{15}.capacity = 30;
    Rooms{15}.Free = true;
    % G.0.14	&	30	
    Rooms{16}.designation = 'G.0.14';
    Rooms{16}.capacity = 30;
    Rooms{16}.Free = true;
    % G.0.15	&	30	
    Rooms{17}.designation = 'G.0.15';
    Rooms{17}.capacity = 30;
    Rooms{17}.Free = true;
    % G.0.16	&	50
    Rooms{18}.designation = 'G.0.16';
    Rooms{18}.capacity = 50;
    Rooms{18}.Free = true;
    % G.0.24	&	81	
    Rooms{19}.designation = 'G.0.24';
    Rooms{19}.capacity = 81;
    Rooms{19}.Free = true;
    % G.1.03	&	50	
    Rooms{20}.designation = 'G.1.03';
    Rooms{20}.capacity = 50;
    Rooms{20}.Free = true;
    % G.1.04	&	45
    Rooms{21}.designation = 'G.1.04';
    Rooms{21}.capacity = 45;
    Rooms{21}.Free = true;
    % G.1.13	&	45
    Rooms{22}.designation = 'G.1.13';
    Rooms{22}.capacity = 45;
    Rooms{22}.Free = true;
    % G.1.15	&	79
    Rooms{23}.designation = 'G.1.15';
    Rooms{23}.capacity = 79;
    Rooms{23}.Free = true;
    % G.1.18	&	40
    Rooms{24}.designation = 'G.1.18';
    Rooms{24}.capacity = 40;
    Rooms{24}.Free = true;
    % G.2.06	&	50
    Rooms{25}.designation = 'G.2.06';
    Rooms{25}.capacity = 50;
    Rooms{25}.Free = true;
    % G.2.07	&	50
    Rooms{26}.designation = 'G.2.07';
    Rooms{26}.capacity = 50;
    Rooms{26}.Free = true;
    % G.2.08	&	50
    Rooms{27}.designation = 'G.2.08';
    Rooms{27}.capacity = 50;
    Rooms{27}.Free = true;
    % G.2.09	&	50
    Rooms{28}.designation = 'G.2.09';
    Rooms{28}.capacity = 50;
    Rooms{28}.Free = true;
    % G.2.10	&	45
    Rooms{29}.designation = 'G.2.10';
    Rooms{29}.capacity = 45;
    Rooms{29}.Free = true;
    % G.2.21	&	48
    Rooms{30}.designation = 'G.2.21';
    Rooms{30}.capacity = 48;
    Rooms{30}.Free = true;
end

% Rooms - Cell array containing room info
function printRooms(Rooms)
    totalCapacity = 0;

    for i = 1 : length(Rooms)
        Rooms{i}.designation
        Rooms{i}.capacity
        totalCapacity = totalCapacity + Rooms{i}.capacity;
    end
    totalCapacity
end


% Lotfi & Cerveny approach
function timetable = assignRooms(Data, timetable)
    % Algorithm:
    %   For each examination period, sort the available rooms in descending 
    %   order by capacity; sort the examinations in descending order by 
    %   number of participants. The largest exam is assigned to the smallest 
    %   room with sufficient capacity to hold the students. If no room is
    %   large enough, then the largest room is filled, and remaining students 
    %   are assigned to a different room.
    %///
    
    % Exam's rooms
    ExamRooms = cell(1, length(Data.ExamCounts));
    
    %   For each examination period, sort the available rooms in descending 
    %   order by capacity
    for p = 1 : timetable.NumPeriods
        % Period rooms
        PeriodRooms = timetable.Rooms;
        % Get capacities
        capacities = zeros(1, length(PeriodRooms));
        for i = 1 : length(PeriodRooms)
            capacities(i) = PeriodRooms{i}.capacity;
        end
        [CapRooms, IdxRooms] = sort(capacities, 'descend');
        % sort the examinations in descending order by number of participants. 
        examList = timetable.Periods{p};
        numExams = length(examList);
        examParticipants = zeros(1, numExams);
        for i = 1 : numExams
            % Compute s_i, the number of students registered for exam e_i
            si = Data.ExamCounts(examList(i));
            examParticipants(examList(i)) = si;
        end    
        [ExamParticipants, IdxExams] = sort(examParticipants, 'descend');
        % The largest exam is assigned to the smallest room with sufficient 
        % capacity to hold the students. If no room is large enough, 
        % then the largest room is filled, and remaining students 
        % are assigned to a different room.
        largestExamIdx = 1;
        remainder = 0;
        while (largestExamIdx <= numExams) 
            capacity = ExamParticipants(largestExamIdx);
            remainder = capacity;
            
            exam = IdxExams(largestExamIdx);
            %Data.Classes{exam}
            
            while (true) 
                [smallestCapRoomIdx, bigEnough, PeriodRooms, assignedRoomCapacity] = ...
                    allocateSmallestRoom(PeriodRooms, IdxRooms, CapRooms, remainder, Data, exam);
                % Assign room
                ExamRooms{IdxExams(largestExamIdx)} = ...
                    [ExamRooms{IdxExams(largestExamIdx)} '-' PeriodRooms{IdxRooms(smallestCapRoomIdx)}.designation]; 
                
                %ExamRooms{IdxExams(largestExamIdx)}
                
                %pause
                
                if (bigEnough)
                    break;
                end
                remainder = remainder - assignedRoomCapacity;
            end
            largestExamIdx = largestExamIdx + 1;
        end
    end
    % Print ExamRooms
    fprintf('Print exam rooms\n');
    for i = 1 : length(ExamRooms)
        fprintf('Exam: ');
        Data.Classes{i}
        fprintf('\nRooms: ');
        ExamRooms{i}
    end
    %pause
end

% Get the smallest room with sufficient capacity to hold the students
function [smallestCapRoomIdx, bigEnough, Rooms, assignedRoomCapacity] = ...
    allocateSmallestRoom(Rooms, IdxRooms, CapRooms, capacity, Data, exam)
    
    for i = length(CapRooms) : -1 : 1
        if (Rooms{IdxRooms(i)}.Free == true && CapRooms(i) >= capacity)
            smallestCapRoomIdx = i;
            bigEnough = true;
            % Mark room as not free
            Rooms{IdxRooms(i)}.Free = false;
            assignedRoomCapacity = CapRooms(i);
            return;
        end
    end
    % If it reach here, there's no room with suficient capacity
    fprintf('There''s no room with suficient capacity for exam %s - capacity = %d !!\n', ...
        Data.Classes{exam}, capacity);
    %pause;
    % Return the greatest free room
    for i = 1 : length(CapRooms)
        if (Rooms{IdxRooms(i)}.Free == true)
            smallestCapRoomIdx = i;
            bigEnough = false;
            % Mark room as not free
            Rooms{IdxRooms(i)}.Free = false;
            assignedRoomCapacity = CapRooms(i);
            return;
        end
    end
    fprintf('There''s not sufficient rooms!!!');
end

%////////////////////////////////////////////////////////////////////////

% Verify Capacity constraint for a timetable period.
% sum {i=1}^{|E|}  a_ip s_i <= S, for all p in P 
% s_i is the number of students registered for exam e_i
function capacity = getPeriodCapacityConstraint(Data, timetable, period) 
    capacity = 0;
    % Get exam list
    examList = timetable.Periods{period};
    numExams = length(examList);
    for i = 1 : numExams
        % Compute s_i, the number of students registered for exam e_i
        si = Data.ExamCounts(examList(i));
        capacity = capacity + si; 
    end
end



%////////////////////////////////////////////////////////////////////////
%
% Hybridization with Stochastic Hill Climber
%
%////////////////////////////////////////////////////////////////////////

% Stochastic Hill Climber used in Hybrid GA
function sol = HCAlgorithm(sol, Data)
    % Check of feasible periods
    tmax = 5; 
%     tmax = 2; 
    T = 0.0001;
    initialSolution = sol;
    % SHC algorithm
    Res = SHC(tmax, T, Data, initialSolution, @getRandomNeigh, @evalFunction);
    sol = Res.s;
end

function c = evalFunction(s, data)
    c = f1aux(s); % Evaluate number of clashes
end

function randomNeigh = getRandomNeigh(Data, timetable)
    randomNeigh = [];
    % In order to identify the most promising moves,
    % a clash list, like the one used in the period expansion operator,
    % is maintained. 
    [clashList, periodList] = buildClashList(Data, timetable);
   
    numClashExams = length(clashList);
    randExamsClashList = randi([1, numClashExams], 1, numClashExams);
    feasibleNeigh = false;
    i = 1;
    while (i <= numClashExams && ~feasibleNeigh)
        % Clashed exam and respective period
        clashedExam = clashList(randExamsClashList(i));
        clashPeriod = periodList(randExamsClashList(i));
        % Determine feasible periods for clashed exam i
        feasiblePeriods = getFeasiblePeriods(Data, timetable, clashedExam, clashPeriod);
        numFeasiblePeriods = length(feasiblePeriods);
        if (numFeasiblePeriods ~= 0)
            % Evaluate moves
            % The move which leads to the greatest decrease in the number
            % of clashes is selected and the exam is removed from the
            % clash list. 
            movesCell = cell(1, numFeasiblePeriods);
            idxMin = -1;
            for p = 1 : numFeasiblePeriods
                % compute move cost
                periodToMove = feasiblePeriods(p);
                movesCell{p} = computeMoveClashes(Data, timetable, clashPeriod, ...
                                        clashedExam, periodToMove);
                % Determine the move which produces the lowest number of clashes
                if (idxMin == -1 || movesCell{p}.numClashes < movesCell{idxMin}.numClashes )
                    idxMin = p;
                end
            end
            % If the min move reduces the number of clashes (of the actual
            % solution), then accept this move.
            randomNeigh = movesCell{idxMin}.timetable;
            t = computeNumClashes(Data, randomNeigh);
            if (randomNeigh.NumClashes ~= t.NumClashes)
                disp('Error')
                t.NumClashes
                randomNeigh.NumClashes
                pause
            end
            
            if (~verifyTimetable(Data, randomNeigh))
                fprintf('randomNeigh is not valid\n');
                pause
            else
                return;
            end
        end
        i = i+1;
    end
    disp('Can''t find feasible random neighbor')
%     pause 
end

function feasiblePeriods = getFeasiblePeriods(Data, timetable, ...
    clashedExam, clashPeriod)
    feasiblePeriods = [];
    % Determine possible periods for clashed exam
    for period = 1 : timetable.NumPeriods
        if (period ~= clashPeriod && isFeasibleExamPeriod(Data, timetable, period, clashedExam))
            feasiblePeriods = [feasiblePeriods period];
        end
    end
end

function [clashList, periodList] = buildClashList(Data, timetable)
    clashList = [];
    periodList = [];
    examListP1 = timetable.Periods{1};
    for p = 2 : timetable.NumPeriods
        examListP2 = timetable.Periods{p};
        if (~(p == 7 || p == 13 || p == 19))
            numExamsP1 = length(examListP1);
            numExamsP2 = length(examListP2);
            marked1 = zeros(1, length(examListP1));
            marked2 = zeros(1, length(examListP2));
            for i = 1 : numExamsP1
                for j = 1 : numExamsP2
                    numStudents = Data.ConflictMatrix(examListP1(i), examListP2(j));
                    if (numStudents > 0) % There's a clash between exams
                        marked1(i) = 1;
                        marked2(j) = 1;
                    end
                end
            end
            % If marked exams aren't yet in the clash list, then copy them.
            [clashList, periodList] = copyExamsClashList(p-1, ...
                periodList, clashList, examListP1(find(marked1 == 1)));
            [clashList, periodList] = copyExamsClashList(p, ...
                periodList, clashList, examListP2(find(marked2 == 1)));
        end
        examListP1 = examListP2;
    end
end

function [clashList periodList] = copyExamsClashList(period, periodList, ...
                                                        clashList, examList)
    if (~isempty(examList))
        for i = 1 : length(examList)
            if (isempty(find(clashList == examList(i))))
                clashList = [clashList examList(i)];
                periodList = [periodList period];
            end
        end
    end
end

function moveInfo = computeMoveClashes(Data, timetable, clashPeriod, ...
                                    clashedExam, period)
    moveInfo.timetable = [];
    moveInfo.numClashes = -1;
    % Verify if clashedExam can move to period (if it's feasible)
    % If it can move, then compute the resulting timetable and number of
    % clashes.
    if (isFeasibleExamPeriod(Data, timetable, period, clashedExam)) 
        % Add clashed exam to period
        timetable.Periods{period} = [timetable.Periods{period} clashedExam];
        % Remove clashed exam from original period
        timetable.Periods{clashPeriod} = remove(clashedExam, timetable.Periods{clashPeriod});
        % Compute number of clashes of new timetable
        moveInfo.timetable = computeNumClashes(Data, timetable);
        moveInfo.numClashes = moveInfo.timetable.NumClashes;
    end
end
%////////////////////////////////////////////////////////////////////////


%////////////////////////////////////////////////////////////////////////
%
% Auxiliary functions
%
%////////////////////////////////////////////////////////////////////////
function Data = removeExamsFromData(Data, ExamsToRemove)
    I = [];
    for i = 1 : length(Data.Classes)
       if (isInExamList(Data.Classes{i}, ExamsToRemove))
           I = [I i];
           % Remove class
           Data.Classes{i} = [];
       end
    end
    % Remove empty cells
    Data.Classes = Data.Classes(~cellfun('isempty', Data.Classes));
    Data.ConflictMatrix(I,:) = [];
    Data.ConflictMatrix(:,I) = [];
end

function b = isInExamList(Class, ExamsToRemove)
    b = false;
    for i = 1 : length(ExamsToRemove)
        if (strcmp(ExamsToRemove{i}, Class))
            b = true;
            return;
        end
    end
end

%////////////////////////////////////////////////////////////////////////
%
% Compute number of students for each course given conflict matrix
%
%////////////////////////////////////////////////////////////////////////
function ExamCounts = computeNumberStudents(C)
    % Compute an approximation of total number of students for each exam.
    % The approximation is givem by the maximum number of students in 
    % the conflict matrix row representing conflicts for that exam.
    N = length(C);
    Exams = zeros(1, N);
    for i = 1 : N
        ExamCounts(i) = max(C(i,:)); 
    end
end

%////////////////////////////////////////////////////////////////////////
%
% Objective functions
%
%////////////////////////////////////////////////////////////////////////
% f1(x) = Number of clashes of solution x
% Sum_{i=1}^{|E|-1} Sum_{j=i+1}^{|E|} Sum_{p=1}^{|P|-1} aip aj(p+1) cij
function y = f1aux(x)
    y = x.NumClashes;
end

% f1(x) = Number of clashes of solution x
% Sum_{i=1}^{|E|-1} Sum_{j=i+1}^{|E|} Sum_{p=1}^{|P|-1} aip aj(p+1) cij
function y = f1(x)
    alpha = 1000;
    y = x.NumClashes;
    Min = x.Range(1);
    Max = x.Range(2);
    if (x.NumPeriods > Max)
        y = y + (x.NumPeriods-Max)*alpha;
    elseif (x.NumPeriods < Min)
        y = y + (Min-x.NumPeriods)*alpha;
    end
end


% f2(x) = Number of periods of solution x
function y = f2(x)
    alpha = 10;
    y = x.NumPeriods;
    Min = x.Range(1);
    Max = x.Range(2);
    if (x.NumPeriods > Max)
        y = y + (x.NumPeriods-Max)*alpha;
    elseif (x.NumPeriods < Min)
        y = y + (Min-x.NumPeriods)*alpha;
    end
end

%////////////////////////////////////////////////////////////////////////

%////////////////////////////////////////////////////////////////////////
%
% Print utilities
%
%////////////////////////////////////////////////////////////////////////
function printTimetable(Data, timetable)
    disp('Timetable info');
    for j = 1 : timetable.NumPeriods
        fprintf('\nPeriod %d', j);
        examList = timetable.Periods{j};
        Data.Classes(examList)
    end
end

function printFinalPopulation(Data, Pop)
    disp('Final population:')
    for i = 1 : length(Pop)
        disp(sprintf('Solution #%d', i));
        x = Pop{i};
        fprintf('# periods = %d\n', x.NumPeriods);
        fprintf('# Clashes = %d\n', x.NumClashes);
        numExams = 0;
        for j = 1 : x.NumPeriods
            fprintf('\nPeriod %d', j);
            examList = x.Periods{j};
            Data.Classes(examList)
            numExams = numExams + length(examList);
        end
        fprintf('# exams = %d\n', numExams);
    end
end

function printPopulation(Data, Pop)
    for i = 1 : Data.NIND
        disp(sprintf('Solution #%d', i));
        x = Pop{i};
        fprintf('# periods = %d\n', x.NumPeriods);
        fprintf('# Clashes = %d\n', x.NumClashes);
        numExams = 0;
        for j = 1 : x.NumPeriods
            examList = x.Periods{j};
            numExams = numExams + length(examList);
        end
    end
end

function numExams = getNumExams(x)
    numExams = 0;
    for j = 1 : x.NumPeriods
        examList = x.Periods{j};
        numExams = numExams + length(examList);
    end
end


%////////////////////////////////////////////////////////////////////////
%
% 1. Initialization
%
%////////////////////////////////////////////////////////////////////////

% 1. Initialization: The population initialization process is
% similar to the reinsertion process of the mutation operator.
% For each chromosome, a timetable
% with a random number of empty periods within the desired
% range is created. Exams are then inserted into randomly selected
% periods in the order determined by the graph coloring
% heuristic, depending on the version of the MOEA. Like
% the mutation operator, when it is not possible to schedule an
% exam without violating any of the hard constraints, a new
% period will be created at the end of the timetable to accommodate
% the exam.
function Pop = initialPopulation(Data)
    % For each chromosome, a timetable with a random number of empty 
    % periods within the desired range is created. 
    Pop = cell(1, Data.NIND);
    for i = 1 : Data.NIND
        Pop{i} = createChromosome(Data);
    end
    
    % Exams are then inserted into randomly selected periods in the order
    % determined by the graph coloring heuristic, depending on the version
    % of the MOEA. Like the mutation operator, when it is not possible to 
    % schedule an exam without violating any of the hard constraints, 
    % a new period will be created at the end of the timetable to 
    % accommodate the exam.
    %
    % Number of exams
    numExams = length(Data.Classes);
    % Generate a random feasible timetable for each element of population
    for i = 1 : Data.NIND
        timetable = Pop{i}; % Timetable contained in chromossome
        fprintf('Generating timetable #%d\n', i);
        UnscheduledExams = 1:numExams;  
        for j = 1 : numExams
            Res = SDheuristicSorting(Data, timetable, UnscheduledExams);
            % Res = SDheuristicSorting1(Data, timetable, UnscheduledExams);
            % Next exam to schedule
            exam = Res.ExamToSchedule;
            if (isempty(Res.ValidPeriods))
                freePeriods = [];
            else
                freePeriods = Res.ValidPeriods(1);
            end
            % Remove exam from UnscheduledExams list
            UnscheduledExams = remove(exam, UnscheduledExams);
            if (length(freePeriods) == 0)
                % It was not possible to schedule this exam without 
                % violating any of the hard constraints. So a new
                % period will be created at the end of the timetable 
                % to accommodate the exam.
                timetable = addPeriodToTimetable(timetable);
                timetable.Periods{timetable.NumPeriods} = [exam];
            else
                % Select a random period from the available periods with free
                % colisions
                idx = randi([1, length(freePeriods)]);
                period = freePeriods(idx);
                if (isFeasibleExamPeriod(Data, timetable, period, exam)) % Just to check
                    timetable.Periods{period} = [timetable.Periods{period} exam];
                else
                    disp('We have a problem!!!')
                    pause;
                end
            end
        end
        if (~verifyTimetable(Data, timetable))
            fprintf('Timetable is not valid\n');
            printTimetable(Data, timetable);
            pause
        end
        timetable = computeNumClashes(Data, timetable);
        Pop{i} = timetable;
    end 
end


function b = verifyTimetable(Data, timetable)
    b = true;
    numExams = getNumExams(timetable);
    if (numExams ~= 80)
        b = false;
        fprintf('Num exams is not 80, but %d\n', numExams);
        b = false;
        return;
    end
    % Verify unicity of exams
    unicityList = [];
    for j = 1 : timetable.NumPeriods
        examList = timetable.Periods{j};
        for k = 1 : length(examList)
            if (~isempty(find(unicityList == examList(k))))
                b = false;
                exam = Data.Classes(examList(k));
                fprintf('%s exam is repeated\n', char(exam));
            else
                unicityList = [unicityList examList(k)];
            end
        end
    end
end

function x = createChromosome(Data)
    % A chromosome encodes a complete and feasible timetable.
    % A period contains a list of exams to be held in the period.
    % For each chromosome, a timetable with a random number of empty 
    % periods within the desired range is created. 
    x.NumPeriods = randi(Data.Range);  % Number of periods 
    x.Range = Data.Range;
    x.Periods = cell(1, x.NumPeriods); % Periods
    x.NumClashes = 0;
    x.Rooms = Data.Rooms;
end

% Compute the number of clashes of solution x
% Sum_{i=1}^{|E|-1} Sum_{j=i+1}^{|E|} Sum_{p=1}^{|P|-1} aip aj(p+1) cij
function x = computeNumClashes(Data, x)
    x.NumClashes = 0;
    examListP1 = x.Periods{1};
    for p = 2 : x.NumPeriods
        examListP2 = x.Periods{p};
        if (~(p == 7 || p == 13 || p == 19))
            numExamsP1 = length(examListP1);
            numExamsP2 = length(examListP2);
            for i = 1 : numExamsP1
                for j = 1 : numExamsP2
                    numStudents = Data.ConflictMatrix(examListP1(i), examListP2(j));
                    if (numStudents > 0)
                        x.NumClashes = x.NumClashes + numStudents;
                    end
                end
            end
        end
        examListP1 = examListP2;
    end
end
 
function numClashes = getNumberClashesPeriod(Data, timetable, period)
    numClashes = 0;
    if (period == 6 || period == 12 || period == 18 || period == 24)
        return;
    end
    examListP1 = timetable.Periods{period};
    examListP2 = timetable.Periods{period+1};
    numExamsP1 = length(examListP1);
    numExamsP2 = length(examListP2);
    for i = 1 : numExamsP1
        for j = 1 : numExamsP2
            numStudents = Data.ConflictMatrix(examListP1(i), examListP2(j));
            if (numStudents > 0)
                numClashes = numClashes + numStudents;
            end
        end
    end
end

% Saturation Degree (SD) graph coloring heuristic: Exams with the 
% fewest valid periods, in terms of satisfying the hard constraints, remaining
% in the timetable are reinserted first.
function Res = SDheuristicSorting(Data, timetable, UnscheduledExams)  
    numExams = length(UnscheduledExams);
    validPeriods = cell(1, numExams);
    for i = 1 : numExams
        % Determine possible periods for exam i
        exam = UnscheduledExams(i);
        P = [];
        for period = 1 : timetable.NumPeriods
             if (isFeasibleExamPeriod(Data, timetable, period, exam))
                numClashes = computeNumClashesExamPeriod(Data, timetable, period, exam);
                %--
                % Extended SD
                %--
                if (numClashes <= 70) 
                    P = [P period];
                end
            end
        end
        validPeriods{i} = P;
    end
    ValidPeriodsCounts = zeros(1, length(validPeriods));
    for i = 1 : length(validPeriods)
        ValidPeriodsCounts(i) = length(validPeriods{i});
    end
    % Sort exams with descending order of fewest valid periods
    [V, I] = sort(ValidPeriodsCounts);
    
    if (V(1) == 0)
        disp('Didn''t found any feasible period');
        idx = randi([1, length(UnscheduledExams)]);
        Res.ExamToSchedule = UnscheduledExams(idx);
        Res.ValidPeriods = [];
        return;
    else
        SortedExams = UnscheduledExams(I); 
        numExams = length(SortedExams);
        ValidPeriods = cell(1, numExams);
        for i = 1 : numExams
            ValidPeriods{i} = validPeriods{I(i)};
        end
        % If there's more than one with the same number of valid periods,
        % choose one randomly
        Ir = find(V == V(1));

        idx = randi([1, length(Ir)]);

        % swap exams
        aux = SortedExams(1);
        per = ValidPeriods{1};

        SortedExams(1) = SortedExams(Ir(idx));
        ValidPeriods{1} = ValidPeriods{Ir(idx)};

        SortedExams(Ir(idx)) = aux;
        ValidPeriods{Ir(idx)} = per;

        Res.ExamToSchedule = SortedExams(1);
        Res.ValidPeriods = ValidPeriods{1};
    end
end


% Saturation Degree (SD) graph coloring heuristic: Exams with the 
% fewest valid periods, in terms of satisfying the hard constraints, remaining
% in the timetable are reinserted first.
function Res = SDheuristicSorting1(Data, timetable, UnscheduledExams)  
    numExams = length(UnscheduledExams);
    validPeriods = cell(1, numExams);
    for i = 1 : numExams
        % Determine possible periods for exam i
        exam = UnscheduledExams(i);
        P = [];
        for period = 1 : timetable.NumPeriods
             if (isFeasibleExamPeriod(Data, timetable, period, exam))
                 P = [P period];
            end
        end
        validPeriods{i} = P;
    end
    ValidPeriodsCounts = zeros(1, length(validPeriods));
    for i = 1 : length(validPeriods)
        ValidPeriodsCounts(i) = length(validPeriods{i});
    end
    % Sort exams with descending order of fewest valid periods
    [V, I] = sort(ValidPeriodsCounts);
    
    if (V(1) == 0)
        disp('Didn''t found any feasible period');
        idx = randi([1, length(UnscheduledExams)]);
        Res.ExamToSchedule = UnscheduledExams(idx);
        Res.ValidPeriods = [];
        return;
    else
        SortedExams = UnscheduledExams(I); 
        numExams = length(SortedExams);
        ValidPeriods = cell(1, numExams);
        for i = 1 : numExams
            ValidPeriods{i} = validPeriods{I(i)};
        end
        % If there's more than one with the same number of valid periods,
        % choose one randomly
        Ir = find(V == V(1));

        idx = randi([1, length(Ir)]);

        % swap exams
        aux = SortedExams(1);
        per = ValidPeriods{1};

        SortedExams(1) = SortedExams(Ir(idx));
        ValidPeriods{1} = ValidPeriods{Ir(idx)};

        SortedExams(Ir(idx)) = aux;
        ValidPeriods{Ir(idx)} = per;

        Res.ExamToSchedule = SortedExams(1);
        Res.ValidPeriods = ValidPeriods{1};
    end
end

function feasible = isFeasibleExamPeriod(Data, timetable, period, exam)
    % Verify exam period feasibility.
    % Constraint that no student is to be scheduled
    % to take two exams at any one time:
    %   Sum_{i=1}^{|E|-1} Sum_{j=i+1}^{|E|} Sum_{p=1}^{|P|} aip ajp cij = 0
    examList = timetable.Periods{period};
    feasible = true;
    
    % Verify capacity constraint
    % Get period current
    capacity = getPeriodCapacityConstraint(Data, timetable, period);
    % Test if capacity was excedeed if exam was inserted
    if (capacity + Data.ExamCounts(exam) > Data.S)
        feasible = false;
        return;
    end
    
    for i = 1 : length(examList)
        numStudents = Data.ConflictMatrix(examList(i), exam);
        if (numStudents > 0)
            feasible = false;
            return;
        end
    end
end

function timetable = addPeriodToTimetable(timetable)
    % Old timetable
    old = timetable.Periods;
    timetable.Periods = cell(1, timetable.NumPeriods+1); % Periods
    % Copy
    size = timetable.NumPeriods;
    for i = 1 : size
        timetable.Periods{i} = old{i};
    end
    timetable.NumPeriods = timetable.NumPeriods+1;
end

% Remove element 'e' from vector 'List'
function List = remove(e, List)
    if (isempty(e))
        return;
    end
    I = find(List == e);
    List = [List(1:I-1) List(I+1:length(List))];
end


%////////////////////////////////////////////////////////////////////////
%
% Crossover
%
%////////////////////////////////////////////////////////////////////////

% Crossing
function Qt = cross(Data, Qt, crossProb)
    probs = rand(Data.NIND/2, 1);
    k = 0;
    for j = 1 : 2 : Data.NIND,
        k = k+1;
        if (probs(k) <= crossProb)
            timetable1 = Qt{j};
            timetable2 = Qt{j+1};
            % Day-exchange crossover
            [timetable1, timetable2] = dayExchangeCrossover(Data, timetable1, timetable2);
            % Actualize chromosomes
            Qt{j} = timetable1;
            Qt{j+1} = timetable2;
        end
    end
end

% Day-exchange crossover
function [timetable1, timetable2] = dayExchangeCrossover(Data, timetable1, timetable2)
    % In day-exchange crossover, only the best days (excluding
    % Saturdays, since exams scheduled on Saturdays are always
    % clash-free) of chromosomes, selected based on the crossover
    % rate, are eligible for exchange. The best day consists of three
    % periods and is the day with the lowest number of clashes per student.
    % 
    % Compute the "Best day" to exchange between chromosomes.
    day1 = bestDay(Data, timetable1);
    day2 = bestDay(Data, timetable2);
    
    if (~verifyTimetable(Data, timetable1))
        fprintf('[Before exchangeDays] Timetable 1 is not valid\n');
        printTimetable(Data, timetable1);
        pause
    end
    if (~verifyTimetable(Data, timetable2))
        fprintf('[Before exchangeDays] Timetable 2 is not valid\n');
        printTimetable(data, timetable2);
        pause
    end
    
    % Exchange best days between the two chromossomes.
    [timetable1, timetable2] = exchangeDays(Data, timetable1, day1, timetable2, day2);
    
    if (~verifyTimetable(Data, timetable1))
        fprintf('[After exchangeDays] Timetable 1 is not valid\n');
        printTimetable(data, timetable1 );
        pause
    end
    if (~verifyTimetable(Data, timetable2))
        fprintf('[After exchangeDays] Timetable 2 is not valid\n');
        printTimetable(data, timetable2);
        pause
    end
    
    % Actualize number of clashes
    timetable1 = computeNumClashes(Data, timetable1);
    timetable2 = computeNumClashes(Data, timetable2);
end


% Compute the "Best day" to exchange between chromosomes.
% The best day consist of three periods (we exclude saturdays because
% saturdays are always clash-free) and is the day with the lowest number
% of clashes per student.
function bday = bestDay(Data, timetable)
    NumPeriods = timetable.NumPeriods;
    clashesPerStudent = zeros(1, NumPeriods);
    % For each period do:
    %
    %   - Get number of clashes in each period divided by 
    %     the number of students in that period.
    numClashes0 = 0;
    
    for p = 1 : NumPeriods
        numStudents = getNumberStudentsPeriod(Data, timetable, p);
        if (p > 1)
            numClashes0 = getNumberClashesPeriod(Data, timetable, p-1); 
        end
        if (p == NumPeriods)
            numClashes = numClashes0;
        else
            numClashes = numClashes0 + getNumberClashesPeriod(Data, timetable, p);              
        end
%         examList = timetable.Periods{p};
%         Data.Classes(examList)
%         clashesPerStudent(p) = numClashes/numStudents; 
        clashesPerStudent(p) = numClashes; 
        
        % VER MELHOR MEDIDA
    
    end
    
    [V, I] = min(clashesPerStudent); % Considering just one period per day
    bday = I(1);
end


% Exchange best days between the two chromossomes.
function [timetable1, timetable2] = exchangeDays(Data, timetable1, day1, timetable2, day2)
    % Exchange days between timetables.
    % The best day and the day after are exchanged -> SEE THIS
    %
    % The newly inserted day takes the place of a randomly chosen day
    % which is pushed to the end of the timetable.
    periodsToInsertTimetable1Copy = cell(1, 1);
    periodsToInsertTimetable1Copy{1} = timetable1.Periods{day1};
    % Insert day2 in timetable1
    idx1 = randi([1, int32(timetable1.NumPeriods)]);
    fromPeriodIdx = idx1;
    periodsToInsert = cell(1, 1);
    periodsToInsert{1} = timetable2.Periods{day2};
    timetable1 = insertDay(Data, timetable1, fromPeriodIdx, periodsToInsert);
    % To ensure the feasibility of chromosomes after the
    % crossover, duplicated exams are deleted. These exams are
    % removed from the original periods, while the newly inserted
    % periods are left intact.
    timetable1 = removeDuplicatedExams(Data, timetable1, fromPeriodIdx, periodsToInsert);

    numExams = 0;
    for j = 1 : timetable1.NumPeriods
        examList = timetable1.Periods{j};
        numExams = numExams + length(examList);
    end
    if (numExams ~= 80)
        pause
    end
    % Insert day1 in timetable2
    idx2 = randi([1, int32(timetable2.NumPeriods)]);
    fromPeriodIdx = idx2;
    periodsToInsert{1} = periodsToInsertTimetable1Copy{1};
    timetable2 = insertDay(Data, timetable2, fromPeriodIdx, periodsToInsert);
    % To ensure the feasibility of chromosomes after the
    % crossover, duplicated exams are deleted. These exams are
    % removed from the original periods, while the newly inserted
    % periods are left intact.
    timetable2 = removeDuplicatedExams(Data, timetable2, fromPeriodIdx, periodsToInsert);
    numExams = 0;
    for j = 1 : timetable2.NumPeriods
        examList = timetable2.Periods{j};
        numExams = numExams + length(examList);
    end

    if (numExams ~= 80)
        pause
    end
end

function showClassNames(Data, periodExamList)
    Data.Classes(periodExamList)
end

function timetable = insertDay(Data, timetable, fromPeriodIdx, periodsToInsert)
    periodsToRemove = cell(1, 1);
    periodsToRemove{1} = timetable.Periods{fromPeriodIdx};
    timetable.Periods{fromPeriodIdx} = periodsToInsert{1};
    % Add two periods to the timetable
    timetable = addPeriodToTimetable(timetable); 
    timetable.Periods{timetable.NumPeriods} = periodsToRemove{1};
end

% To ensure the feasibility of chromosomes after the
% crossover, duplicated exams are deleted. These exams are
% removed from the original periods, while the newly inserted
% periods are left intact.
function timetable = removeDuplicatedExams(Data, timetable, fromPeriodIdx, periodsToInsert)
    % Form an exam list with the inserted periods
    insertedExams = [periodsToInsert{1}];
    numRepeatedExams = 0;
    
    for i = 1 : timetable.NumPeriods
        % Except the periods just inserted
         if (i ~= fromPeriodIdx) 
            % Remove duplicated exams
            examList = timetable.Periods{i};
            numExams = length(examList);
            % Find exams from this list which are repeated and remove them
            exams = [];
            for j = 1 : numExams
                I = find(insertedExams == examList(j));
                if (~isempty(I))
                    exams = [exams examList(j)];
                    numRepeatedExams = numRepeatedExams+1;
                end
            end
            % Remove exams
            if (~isempty(exams))
                for j = 1 : length(exams)
                    examList = remove(exams(j), examList);
                end
                timetable.Periods{i} = examList;
            end
        end
    end
    if (numRepeatedExams ~= length(insertedExams))
        disp('Bug...')
        pause
    end
end

%////////////////////////////////////////////////////////////////////////
%
% Mutation
%
%////////////////////////////////////////////////////////////////////////
function Qt = mutate(Data, Qt, mutProb, reinsertionRate)
    % For each chromosome selected for mutation based on the mutation rate,
    % the operator removes a number of exams, selected based on the 
    % reinsertion rate, from the chromosome. These exams are then 
    % reinserted into randomly selected periods while maintaining feasibility. 
    probs = rand(Data.NIND, 1);
    for j = 1 : Data.NIND,
        if (probs(j) < mutProb)
            timetable = Qt{j};
            % Remove a number of exams, selected based on the 
            % reinsertion rate, from the chromosome.
            [timetable, UnscheduledExams] = removeExamsForReinsertion(Data, ...
                                               timetable, reinsertionRate);

            % These exams are then reinserted into randomly selected 
            % periods while maintaining feasibility.
            timetable = reinsertExams(Data, timetable, UnscheduledExams);
            timetable = removeEmptyPeriods(Data, timetable);            

            if (~verifyTimetable(Data, timetable))
                fprintf('[After mutate] Timetable is not valid\n');
                printTimetable(Data, timetable);
                pause
            end
            timetable = computeNumClashes(Data, timetable);
             % Actualize chromosome
            Qt{j} = timetable;
        end
    end
    
    Qt = packTimetables(Data, Qt);
    
end


function [timetable, UnscheduledExams] = removeExamsForReinsertion(Data, ...
                                                timetable, reinsertionRate)
    numExams = length(Data.Classes);
    UnscheduledExams = [];
    probs = rand(numExams, 1);
    for j = 1 : Data.NIND
        if (probs(j) < reinsertionRate)
            timetable = removeExam(timetable, j);
            UnscheduledExams = [UnscheduledExams j];
        end
    end
end

function timetable = removeExam(timetable, exam)
    for i = 1 : timetable.NumPeriods
        % Remove exam
        examList = timetable.Periods{i};
        I = find(examList == exam);
        % If found, remove exam
        if (~isempty(I))
            examList = remove(examList(I), examList);
            timetable.Periods{i} = examList;
        end
    end
end

% Reinsertion process of the mutation operator.
% Exams are inserted into randomly selected
% periods in the order determined by the graph coloring
% heuristic, depending on the version of the MOEA. 
% When it is not possible to schedule an
% exam without violating any of the hard constraints, a new
% period will be created at the end of the timetable to accommodate
% the exam.
function timetable = reinsertExams(Data, timetable, UnscheduledExams)
    % Exams are then inserted into randomly selected periods in the order
    % determined by the graph coloring heuristic, depending on the version
    % of the MOEA. Like the mutation operator, when it is not possible to 
    % schedule an exam without violating any of the hard constraints, 
    % a new period will be created at the end of the timetable to 
    % accommodate the exam.
    %
    % Number of exams to schedule
    numExams = length(UnscheduledExams);
    for j = 1 : numExams
        Res = SDheuristicSorting1(Data, timetable, UnscheduledExams);
        % Next exam to schedule
        exam = Res.ExamToSchedule;
        if (isempty(Res.ValidPeriods))
            freePeriods = [];
        else
            freePeriods = Res.ValidPeriods(1);
        end
        % Remove exam from UnscheduledExams list
        UnscheduledExams = remove(exam, UnscheduledExams);

        if (length(freePeriods) == 0)
            % It was not possible to schedule this exam without 
            % violating any of the hard constraints. So a new
            % period will be created at the end of the timetable 
            % to accommodate the exam.
            timetable = addPeriodToTimetable(timetable);
            timetable.Periods{timetable.NumPeriods} = [exam];
        else
            % Select a random period from the available periods with free
            % colisions
            idx = randi([1, length(freePeriods)]);
            period = freePeriods(idx);
            if (isFeasibleExamPeriod(Data, timetable, period, exam)) % Just to check
                timetable.Periods{period} = [timetable.Periods{period} exam];
            else
                disp('We have a problem!!!')
                pause;
            end
        end
    end
end

%////////////////////////////////////////////////////////////////////////
%
% Timetable packing
%
%////////////////////////////////////////////////////////////////////////
% Pack timetables
function Pop = packTimetables(Data, Pop)
%     fprintf('\nPacking timetables...\n');
    for i = 1 : length(Pop)
        x = Pop{i};
        if (x.NumPeriods > Data.MaxPeriods)
            % Timetable period packing
%             fprintf('Packing timetable %d\n', i);
            x = periodPacking(Data, x);
        elseif (x.NumPeriods < Data.MinPeriods)
            % Timetable period expansion
%             fprintf('Expanding timetable %d\n', i);
            x = periodExpansion(Data, x);
        end
        Pop{i} = x;
    end
end


% Timetable period expansion
function timetable = periodExpansion(Data, timetable)
    % Period expansion: The operation first adds empty periods
    % to the end of the timetable such that the timetable
    % length is equal to a random number within the desired range.
    % A clash list, consisting of all exams that are involved in at
    % least one clash, is also maintained. An exam is randomly
    % selected from the clash list and the operation searches in a
    % random order for a period which the selected exam can be
    % rescheduled without causing any clashes while maintaining
    % feasibility. The exam remains intact if no such period exists.
    % The operation ends after one cycle through all exams in the
    % clash list.
    
    % Feasible number of periods
    numPeriods = randi([Data.MinPeriods Data.MaxPeriods]);
    numExtraPeriods = numPeriods - timetable.NumPeriods;
    % add empty periods to the end of the timetable such that the timetable
    % length is equal to a random number within the desired range.
    for i = 1 : numExtraPeriods
        timetable = addPeriodToTimetable(timetable);
    end
    % A clash list, consisting of all exams that are involved in at
    % least one clash, is also maintained. 
    [clashList, periodList] = buildClashList(Data, timetable);
    
    while (~isempty(clashList))
        % An exam is randomly
        % selected from the clash list and the operation searches in a
        % random order for a period which the selected exam can be
        % rescheduled without causing any clashes while maintaining
        % feasibility. The exam remains intact if no such period exists.
        % The operation ends after one cycle through all exams in the
        % clash list.

        % Select an exam randomly from the clash list
        examToScheduleIdx = randi([1 length(clashList)]);
        examToSchedule = clashList(examToScheduleIdx);
        period = periodList(examToScheduleIdx);
        
        % Remove exam from clash list
        clashList = remove(examToSchedule, clashList);
        periodList = remove(period, periodList);
        
        % Get available periods (set of periods where the exam can be
        % rescheduled without causing any clashes while maintaining
        % feasibility)
        numPeriods = timetable.NumPeriods;
        AllPeriodsCapacity = zeros(1, numPeriods);
        % Get available periods capacity
        for p = 1 : numPeriods
            if (p ~= period)
                AllPeriodsCapacity(p) = getPeriodCapacity(Data, ...
                    timetable, p, examToSchedule)
            end
        end
        ValidPeriodsCapacityIdxs = find(AllPeriodsCapacity);
        ValidPeriodsCapacity = AllPeriodsCapacity(ValidPeriodsCapacityIdxs);
        if (~isempty(ValidPeriodsCapacity))
            % the operation searches in a
            % random order for a period which the selected exam can be
            % rescheduled without causing any clashes while maintaining
            % feasibility. The exam remains intact if no such period exists.
            % The operation ends after one cycle through all exams in the
            % clash list. 
            idx = randi([1 length(ValidPeriodsCapacity)]);
            % Period
            p = ValidPeriodsCapacityIdxs(idx);
            % Reschedule exam
            timetable.Periods{p} = [timetable.Periods{p} examToSchedule];
            % Remove exam from original period
            timetable.Periods{period} = remove(examToSchedule, ...
                                            timetable.Periods{period});
            
            if (~verifyTimetable(Data, timetable))
                printTimetable(Data, timetable);
                fprintf('[After expansion] Timetable is not valid\n');
                pause
            end
            if (~isFeasibleTimetable(Data, timetable))
                printTimetable(Data, timetable);
                fprintf('[After expansion] Timetable is not feasible\n');
                pause
            end
        else
%             disp('no valid periods')
%             pause
        end
    end
    % Recompute number of clashes
    timetable = computeNumClashes(Data, timetable)
end


% Get period's available capacity
function capacity = getPeriodCapacity(Data, timetable, period, examToSchedule)
    capacity = 0;
    % Verify if period 'p' is a feasible period for current exam 
    if (isFeasibleExamPeriod(Data, timetable, period, examToSchedule))
        % Compute number of clashes between 'exam' and other exams
        % in periods p-1 and p+1
        numClashes = computeNumClashesExamPeriod(Data, timetable, period, examToSchedule);
        if (numClashes == 0)
            capacity = 1;
        end
    end
 end


function feasible = isFeasibleTimetable(Data, timetable)
    feasible = true;
    for p = 1 : timetable.NumPeriods
        examList = timetable.Periods{p};
        for i = 1 : length(examList)
            feasible = isFeasibleExamPeriod(Data, timetable, p, examList(i));
            if (~feasible)
                disp('Not feasible')
                Data.Classes{examList(i)}
                p
                return;
            end
        end
    end
end

    
% Timetable period packing
function timetable = periodPacking(Data, timetable)
    % Period packing: Starting from the period with the
    % smallest number of students, the operation searches in order
    % of available period capacity, starting from the smallest,
    % for a period which can accommodate exams from the former
    % without causing any clashes while maintaining feasibility.
    % The operation stops when it goes one cycle through all periods
    % without rescheduling any exam or when the timetable
    % length is reduced to a random number within the desired
    % range.
    
    timetable = removeEmptyPeriods(Data, timetable);
    numPeriods = timetable.NumPeriods;
    % Get period with the smallest number of students.
    minPeriod = 1;
    min = getNumberStudentsPeriod(Data, timetable, 1);
    for p = 2 : numPeriods
        count = getNumberStudentsPeriod(Data, timetable, p);
        if (count < min)
            min = count;
            minPeriod = p;
        end
    end
    examList = timetable.Periods{minPeriod};
    % The operation searches in order of available period capacity, 
    % starting from the smallest, for a period which can accommodate exams
    % from the former without causing any clashes while maintaining 
    % feasibility. 
    % The operation stops when it goes one cycle through all periods
    % without rescheduling any exam or when the timetable
    % length is reduced to a random number within the desired range.
    AllPeriodsCapacity = zeros(1, numPeriods);
    % Get available periods capacity
    for p = 1 : numPeriods
        if (p ~= minPeriod)
            AllPeriodsCapacity(p) = getPeriodAvailableCapacity(Data, timetable, p, minPeriod);
        end
    end
    ValidPeriodsCapacityIdxs = find(AllPeriodsCapacity);
    ValidPeriodsCapacity = AllPeriodsCapacity(ValidPeriodsCapacityIdxs);
    
    if (~isempty(ValidPeriodsCapacity))
       % Sort periods by descending order of available period capacity
       [V, I] = sort(ValidPeriodsCapacity);
       % Reschedule an exam 
       LeastAvailablePeriodsIdxs = find(ValidPeriodsCapacity == ValidPeriodsCapacity(1));
       idx = randi([1, length(LeastAvailablePeriodsIdxs)]);
       
       p = ValidPeriodsCapacityIdxs(LeastAvailablePeriodsIdxs(idx));
       
       [numExams, ExamList] = getPeriodAvailableCapacity(Data, timetable, ...
           p, minPeriod);
       % Reschedule exams in the returned exam list from period 'minPeriod'
       % and insert them in period 'p'
       examListMinPeriod = timetable.Periods{minPeriod};
       for i = 1 : length(ExamList)
           exam = ExamList(i);
           % Remove exam from minPeriod
           examListMinPeriod = remove(exam, examListMinPeriod);
           % Reschedule exam
           timetable.Periods{p} = [timetable.Periods{p} exam];
       end
       timetable.Periods{minPeriod} = examListMinPeriod;
       
       if (isempty(examListMinPeriod))
           % Remove minPeriod
            timetable = removeEmptyPeriods(Data, timetable);
        end
    end
    % Recompute number of clashes
    timetable = computeNumClashes(Data, timetable);
end

function timetable = removeEmptyPeriods(Data, timetable)
    timetable.Periods = timetable.Periods(~cellfun('isempty', timetable.Periods));
    timetable.NumPeriods = length(timetable.Periods);
end


% Get period's available capacity
function [capacity, exams] = getPeriodAvailableCapacity(Data, timetable, p, minPeriod)
    % p period's capacity is defined as the number of exams from 
    % period 'minPeriod' which can accommodated without causing 
    % any clashes while maintaining feasibility. 
    capacity = 0;
    exams = [];
    examList = timetable.Periods{minPeriod};
    for i = 1 : length(examList)
        exam = examList(i);
        % Verify if period 'p' is a feasible period for current exam 
        if (isFeasibleExamPeriod(Data, timetable, p, exam))
            % Compute number of clashes between 'exam' and other exams
            % in periods p-1 and p+1
            numClashes = computeNumClashesExamPeriod(Data, timetable, p, exam);
            if (numClashes == 0)
                capacity = capacity+1;
                exams = [exams exam];
            end
        end
    end
 end

% Compute number of clashes between 'exam' and other exams
% in periods p-1 and p+1
function numClashes = computeNumClashesExamPeriod(Data, timetable, p, exam)
    if (p == 6 || p == 12 || p == 18)
        numClashes1 = 0;
    else
        numClashes1 = 0;
        if (p+1 <= timetable.NumPeriods)
            examListP1 = timetable.Periods{p+1};
        else
            examListP1 = [];
        end
        numExamsP1 = length(examListP1);
        for i = 1 : numExamsP1
            numStudents = Data.ConflictMatrix(examListP1(i), exam);
            if (numStudents > 0)
                numClashes1 = numClashes1 + numStudents;
            end
        end
    end
    
    if (p == 7 || p == 13 || p == 19)
        numClashes0 = 0;
    else
        numClashes0 = 0;
        if (p-1 >= 1)
            examListP0 = timetable.Periods{p-1};
        else
            examListP0 = [];
        end
        numExamsP0 = length(examListP0);
        for i = 1 : numExamsP0
            numStudents = Data.ConflictMatrix(examListP0(i), exam);
            if (numStudents > 0)
                numClashes0 = numClashes0 + numStudents;
            end
        end
    end
    numClashes = numClashes0 + numClashes1;
end
    

% Get number of students in a period
function studentCount = getNumberStudentsPeriod(Data, timetable, period)
    examList = timetable.Periods{period};
    studentCount = 0;
    for i = 1 : length(examList)
        studentCount = studentCount + Data.ExamCounts(examList(i));
    end
end
