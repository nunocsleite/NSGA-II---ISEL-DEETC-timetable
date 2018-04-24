%-------------------------------------------------------------
%  Final Project of OMBAE PhD course, Winter semester 2011/12
%
%  Author: Nuno Miguel Leite, No. 59734/D
%-------------------------------------------------------------
function generateExcelTimetables(fileName, Data, timetable)
    % timetable parameter is the merged timetable

    filenames = {'LEETC'; 'LEIC'; 'LERCM'; 'MEIC'; 'MEET'};
    
    ExamsLEETC = {'ALGA', 'Pg', 'AM1', 'FAE', 'ACir', ...
                   'POO', 'AM2', 'LSD', 'E1', 'MAT', ...
                   'PE', 'ACp', 'EA', 'E2', 'SS', ...
                   'RCp', 'PICC/CPg', 'PR', 'FT', 'SEAD1', ...
                   'ST', 'RCom', 'SEAD2', 'RI', 'SE1', 'AVE', 'SCDig', 'SOt', ...
                   'PI', 'SCDist', 'EGP', 'OGE', 'SG'};

               
    ExamsLEIC = { 'ALGA', 'Pg', 'LSD', 'M1', 'Elctr', ...
                'POO', 'PE', 'ACp', 'M2' 'PSC', ...
                'AC2', 'AED', 'Com', 'CG', 'SI1', 'LC', 'PF', 'EGP', 'OGE', 'SG', ...
                'RCp', 'AVE', 'SOi', ...
                'PI', 'SI2', 'PC', 'SI', 'RI', 'SE1', 'Cpl', ...
                'SD' };

    ExamsLERCM = { 'ALGA', 'Pg', 'AM1', 'F1', 'ITI', ...
                'POO', 'PE', 'AM2', 'F2', 'PDSr', ...
                'SCDig', 'MNO', 'AIEC', 'PICC/CPg', 'CSDin', ...
                'RCp', 'CSM', 'FIA', 'MSr', 'SOt', ...
                'BD', 'CGAV', 'AA', 'RI', 'SCDist', ...
                'ES', 'EGP', 'OGE', 'RSCM', 'PCM', 'PIV'              
        };
    
    ExamsMEIC = { 'SI', 'SD', 'ES', 'RI', 'SE1', 'Cpl', 'CCD', 'CSO', ...
                                'CSI', 'RSCM', 'CAC', 'CIA', 'AA', ...
        'ASI', 'GSI', 'PSTR', 'IRS', 'IS', 'EGP', 'OGE', 'SG'     };
    
    ExamsMEET = { 'AVE', 'SE1', 'RI', 'SEAD2', 'ST', 'RCom', 'PIV', ...
                        'SCDig', 'SEAS', 'CEE', 'OE', 'PRC', 'RSCM', 'SET', ... 
                        'Ant', 'CMov', 'STBL', 'CSDist', ...
                  'BD', 'CSM', 'PI', 'SCDist', 'EGP', 'OGE', 'SG', 'SCO', ...
                    'PCI', 'IRS', 'RDC', 'SEADI', 'RM', 'STDS/PSTR'
        };
    
    ExamsNames = cell(5,1);
    ExamsNames{1} = ExamsLEETC;
    ExamsNames{2} = ExamsLEIC;               
    ExamsNames{3} = ExamsLERCM;
    ExamsNames{4} = ExamsMEIC;    
    ExamsNames{5} = ExamsMEET;     
    
    mergedTimetable = [];
        
    for i = 1 : length(ExamsNames)
        periods = getClassPeriods(Data, timetable, ExamsNames{i});
        numPeriods = timetable.NumPeriods;
        % Build timetable
        courseTimetable = buildTimetable(Data, periods, ExamsNames{i}, ...
                                                filenames{i}, numPeriods);
        courseTimetable = computeNumClashes(Data, courseTimetable);
        verifyTimetable(Data, courseTimetable)
        
        printTimetableFile(strcat(fileName, '_', filenames{i}, '.csv'), ...
            ExamsNames{i}, periods, numPeriods, courseTimetable.NumClashes);
        
        % Merge timetables
        mergedTimetable = merge(Data, mergedTimetable, courseTimetable);
    end
    
    % For DEBUG purposes
    mergedTimetable.name = 'Merged';
    mergedTimetable = computeNumClashes(Data, mergedTimetable);
    disp('Final Merged timetable')
    printTimetable(Data, mergedTimetable)
    verifyTimetable(Data, mergedTimetable)
    
    disp('Comparing number of clashes')
    timetable.NumClashes
    mergedTimetable.NumClashes
    
    disp('Compute number of clashes again')
    timetable = computeNumClashes(Data, timetable);
    timetable.NumClashes
    mergedTimetable.NumClashes
end

% ////////////////////////////////////////////////////////////////////////
% Merge timetables
function timetable = merge(Data, timetable, toMerge)
    if (isempty(timetable))
        timetable = toMerge;
    else
        % Copy exams
        for p = 1 : timetable.NumPeriods
            examList1 = timetable.Periods{p};
            examList2 = toMerge.Periods{p};
            examList1 = mergeExams(examList1, examList2);
            timetable.Periods{p} = examList1;
        end
    end
end

% Merge exams
function examList1 = mergeExams(examList1, examList2)
    numExams2 = length(examList2);
    for i = 1 : numExams2
        % If it doesn't exist, copy it
        if (isempty(find(examList1 == examList2(i))))
            examList1 = [examList1 examList2(i)];
        end
    end
end
% ////////////////////////////////////////////////////////////////////////

% ////////////////////////////////////////////////////////////////////////
function b = verifyTimetable(Data, timetable)
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
    feasible = isFeasibleTimetable(Data, timetable)
end

function feasible = isFeasibleTimetable(Data, timetable)
    feasible = true;
    for p = 1 : timetable.NumPeriods
        examList = timetable.Periods{p};
        for i = 1 : length(examList)
            feasible = isFeasibleExamPeriod(Data, timetable, p, examList(i));
            if (~feasible)
                disp('Not feasible')
                pause
            end
        end
    end
end

function feasible = isFeasibleExamPeriod(Data, timetable, period, exam)
    % Verify exam period feasibility.
    % Constraint that no student is to be scheduled
    % to take two exams at any one time:
    %   Sum_{i=1}^{|E|-1} Sum_{j=i+1}^{|E|} Sum_{p=1}^{|P|} aip ajp cij = 0
    examList = timetable.Periods{period};
    feasible = true;
    InfeasibleExamPairs = [];
    numInfeasibleExamPairs = 0;
    for i = 1 : length(examList)
        numStudents = Data.ConflictMatrix(examList(i), exam);
        if (numStudents > 0)
            feasible = false;
                numInfeasibleExamPairs = numInfeasibleExamPairs+1;
                InfeasibleExamPairs = [InfeasibleExamPairs; ...
                    [Data.Classes(examList(i)) Data.Classes(exam)] ];
                disp(sprintf('exam %d - Number of infeasible exam pairs: %d\n\n', exam, numInfeasibleExamPairs))
                InfeasibleExamPairs
                pause
            return;
        end
    end
end
% ////////////////////////////////////////////////////////////////////////

function timetable = buildTimetable(Data, periods, ExamsNames, name, numPeriods)
    % Create cromosome
    timetable = createChromosome(Data, name, numPeriods);
    for p = 1 : length(periods)
        examName = ExamsNames{p};
        examIdx = getExamIndex(Data, examName);
        if (periods(p) ~= 0)
            timetable.Periods{periods(p)} = [timetable.Periods{periods(p)} examIdx];
        end
    end
end

function printTimetable(Data, timetable)
    disp('Timetable info');
    fprintf('timetable: %s\n', timetable.name);
    fprintf('# periods = %d\n', timetable.NumPeriods);
    fprintf('# Clashes = %d\n', timetable.NumClashes);
    numExams = 0;
    for j = 1 : timetable.NumPeriods
        fprintf('\nPeriod %d', j);
        examList = timetable.Periods{j};
        Data.Classes(examList)
        numExams = numExams + length(examList);
    end
    fprintf('# exams = %d\n\n', numExams);
end


%////////////////////////////////////////////////////////////////////////
%
% Compute the number of clashes of solution x
% Sum_{i=1}^{|E|-1} Sum_{j=i+1}^{|E|} Sum_{p=1}^{|P|-1} aip aj(p+1) cij
%
%////////////////////////////////////////////////////////////////////////
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

%////////////////////////////////////////////////////////////////////////
%
% Chromosome representation
%
%////////////////////////////////////////////////////////////////////////
function x = createChromosome(Data, name, numPeriods)
    % A chromosome encodes a complete and feasible timetable.
    % A period contains a list of exams to be held in the period.
    % For each chromosome, a timetable with a random number of empty 
    % periods within the desired range is created. 
    x.NumPeriods = numPeriods;  % Number of periods 
    x.Periods = cell(1, x.NumPeriods); % Periods
    x.NumClashes = 0;
    x.name = name;
end

function periods = getClassPeriods(Data, timetable, Classes)
    periods = zeros(length(Classes),1);
    for i = 1 : length(periods)
        class = Classes{i};
        classIdx = getExamIndex(Data, class);
        per = 0;
        % Get period number for this class
        for p = 1 : timetable.NumPeriods
            examList = timetable.Periods{p};
            for e = 1 : length(examList)
                if (examList(e) == classIdx)
                    per = p;
                    break;
                end
            end
        end
        periods(i) = per;
    end
end

function examIdx = getExamIndex(Data, examName)
    examIdx = -1; % Error, not found
    for i = 1 : length(Data.Classes)
        if (strcmp(Data.Classes{i}, examName))
            examIdx = i;
            return;
        end
    end
end

% Print timetable in tabular form to the console.
function printTimetableFile(filename, ExamsLEETC, Values, P, numClashes)
    fid = fopen(filename, 'wt');
    fprintf(fid, 'Num. Clashes = %d\n\n', numClashes);
    numClasses = length(ExamsLEETC);
    numPeriods = P;
    printPeriodsRow(fid, numPeriods)
    for i = 1 : numClasses
        printRow(fid, i, ExamsLEETC, Values, numPeriods);
    end
    fclose(fid);
end

function printPeriodsRow(fid, numPeriods)
    for i = 1 : numPeriods
        fprintf(fid, ';%d', i);
    end
    fprintf(fid, '\n');
end

function printRow(fid, i, ExamsLEETC, Values, numPeriods)
    className = ExamsLEETC{i};
    fprintf(fid, '%s;', className);
    for p = 1 : numPeriods
        if (Values(i) == p)
            c = 'x;';
        else
            c = ';';
        end
        fprintf(fid, '%s', c);
    end
    fprintf(fid, '\n');
end
