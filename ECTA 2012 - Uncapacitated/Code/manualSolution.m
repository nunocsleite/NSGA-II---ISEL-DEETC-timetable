%-------------------------------------------------------------
%  Final Project of OMBAE PhD course, Winter semester 2011/12
%
%  Author: Nuno Miguel Leite, No. 59734/D
%-------------------------------------------------------------
function manualSolution()
    % ISEL/DEETC Examination Timetabling Problem
    % 
    clc; % Clear screen
     
    fileToRead = 'ManualSolution.xls';
    %
    % ISEL/DEETC department test data
    %
    Data = importfile('Relacoes UC 0910i.xls');

    Data.ExamsToRemove = {'LIC', 'SEAD2', 'AC2', 'FIA', 'PCM', 'CIA', 'PSTR', ...
        'SEAS', 'CEE', 'PRC', 'SEADI', 'RM', 'STDS'}; % PRC does not have exam
    Data = removeExamsFromData(Data, Data.ExamsToRemove);
    
    ExamCounts = computeNumberStudents(Data.ConflictMatrix);
    Data.ExamCounts = ExamCounts;
     % Periods range definition. There are 3 period in each day
    % and one period at saturday.
    NumberWeeks = 3;
    % Number of periods.
    Data.NumPeriods = NumberWeeks*6; % 1 period per day
    
    sheetName = 'LEETC';
    % Import the file
    [timetableLEETC] = importSheetData(Data, fileToRead, sheetName);

    sheetName = 'LEIC';
    % Import the file
    [timetableLEIC] = importSheetData(Data, fileToRead, sheetName);
    
    sheetName = 'LERCM';
    % Import the file
    [timetableLERCM] = importSheetData(Data, fileToRead, sheetName);
    
    sheetName = 'MEIC';
    % Import the file
    [timetableMEIC] = importSheetData(Data, fileToRead, sheetName);
    
    sheetName = 'MEET';
    % Import the file
    [timetableMEET] = importSheetData(Data, fileToRead, sheetName);
    
    timetables = cell(1, 5);
    timetables{1} = timetableLEETC;
    timetables{2} = timetableLEIC;
    timetables{3} = timetableLERCM;
    timetables{4} = timetableMEIC;
    timetables{5} = timetableMEET;
    timetable = [];
    fileNames = {'LEETC.csv'; 'LEIC.csv'; 'LERCM.csv'; 'MEIC.csv'; 'MEET.csv'};
    
    for i = 1 : 5
        fileNames{i}
        printTimetableFile(fileNames{i}, Data, timetables{i});
        disp('Done')
        pause
        timetables{i} = computeNumClashes(Data, timetables{i});
        timetables{i}.NumClashes
        printTimetable(Data, timetables{i})
        verifyTimetable(Data, timetables{i})
        % Merge timetables
        timetable = merge(Data, timetable, timetables{i});
    end
    timetable.name = 'Merged';
    timetable = computeNumClashes(Data, timetable);
    printTimetable(Data, timetable)
    verifyTimetable(Data, timetable)
end

% Print timetable in tabular form to the console.
function printTimetableFile(fileName, Data, timetable)
    fid = fopen(fileName, 'wt');
    numClasses = length(timetable.Classes);
    numPeriods = timetable.NumPeriods;
    printPeriodsRow(fid, numPeriods)
    for i = 1 : numClasses
        printRow(fid, Data, i, timetable);
    end
    fclose(fid);
end

function printPeriodsRow(fid, numPeriods)
    for i = 1 : numPeriods
        fprintf(fid, ';%d', i);
    end
    fprintf(fid, '\n');
end

function printRow(fid, Data, i, timetable)
    numPeriods = timetable.NumPeriods;
    className = timetable.Classes{i};
    examIdx = getExamIndex(Data, className);
    fprintf(fid, '%s;', className);
    for p = 1 : numPeriods
        examList = timetable.Periods{p};
        if (~isempty(find(examList == examIdx)))
            c = 'x;';
        else
            c = ';';
        end
        fprintf(fid, '%s', c);
    end
    fprintf(fid, '\n');
end

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

function [timetable] = importSheetData(Data, fileToRead, sheetName)
    [numbers, strings] = xlsread(fileToRead, sheetName);
    if ~isempty(numbers)
        newData.data =  numbers;
    end
    if ~isempty(strings)
        newData.textdata =  strings;
    end

    % Create cromosome
    timetable = createChromosome(Data, sheetName);
    % Fill with data
    classes = newData.textdata(:,1);
    table = newData.textdata(:, 2:19); % 1st epoch
    timetable.Classes = classes;
    
    numExams = length(Data.Classes);
    scheduled = zeros(1, numExams);
    
    for p = 1 : timetable.NumPeriods
        examsIdxs = findExamIdxs(table(:, p));
        for i = 1 : length(examsIdxs)
            examName = classes(examsIdxs(i));
            % If not removed and not scheduled yet (there are repeated markings 
            % in the manual timetable), then schedule it
           if (~isInExamList(examName, Data.ExamsToRemove))
               examIdx = getExamIndex(Data, examName);
               if (scheduled(examIdx) == 0)
                   scheduled(examIdx) = 1;
                   timetable.Periods{p} = [timetable.Periods{p} examIdx];
               end
           end
        end
    end
end


%////////////////////////////////////////////////////////////////////////
%
% Auxiliary functions
%
%////////////////////////////////////////////////////////////////////////
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

function examsIdxs = findExamIdxs(tablePeriod)
    examsIdxs = [];
    for i = 1 : length(tablePeriod)
        if (strcmp(tablePeriod(i), 'x'))
            examsIdxs = [examsIdxs i];
        end
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

function numExams = getNumExams(x)
    numExams = 0;
    for j = 1 : x.NumPeriods
        examList = x.Periods{j};
        numExams = numExams + length(examList);
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
% Chromosome representation
%
%////////////////////////////////////////////////////////////////////////
function x = createChromosome(Data, name)
    % A chromosome encodes a complete and feasible timetable.
    % A period contains a list of exams to be held in the period.
    % For each chromosome, a timetable with a random number of empty 
    % periods within the desired range is created. 
    x.NumPeriods = Data.NumPeriods;  % Number of periods 
    x.Periods = cell(1, x.NumPeriods); % Periods
    x.NumClashes = 0;
    x.name = name;
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


