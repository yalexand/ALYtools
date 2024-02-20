function cell2csv(fileName, cellArray, numericPrecision, delimiter, filePermissions, excelYear, decimal)
% % Writes cell array content into a *.csv file.
% % 
% % CELL2CSV(fileName, cellArray[, delimiter, excelYear, decimal])
% %
% % fileName     = Name of the file to save. [ e.g. 'text.csv' ]
% % cellArray    = Name of the Cell Array where the data is in
% % 
% % optional:
% % numericPrecision = determines the output precision level of numbers using standard precision syntax
% %                    (eg. '%.12f')
% % delimiter    = sign separating the values (default = ',')
% % filePermissions = opens the file fileName in the mode specified by
% %     filePermissions:
% %     'r'     open file for reading
% %     'w'     open file for writing; discard existing contents
% %     'a'     open or create file for writing; append data to end of file
% %     'r+'    open (do not create) file for reading and writing
% %     'w+'    open or create file for reading and writing; discard existing
% %             contents
% %     'a+'    open or create file for reading and writing; append data to end of
% %             file             
% %     'W'     open file for writing without automatic flushing
% %     'A'     open file for appending without automatic flushing
% % excelYear    = depending on the Excel version, the cells are put into
% %                quotes before they are written to the file. The delimiter
% %                is set to semicolon (;)  (default = 1997 which does not change delimiter to semicolon ;)
% % decimal      = defines the decimal delimiter (default = '.')
% %
% %         by Sylvain Fiedler, KA, 2004
% % updated by Sylvain Fiedler, Metz, 06
% % fixed the logical-bug, Kaiserslautern, 06/2008, S.Fiedler
% % added the choice of decimal delimiter, 11/2010, S.Fiedler
% % modfiedy and optimized by Jerry Zhu, June, 2014, jerryzhujian9@gmail.com
% % now works with empty cells, numeric, char, string, row vector, and logical cells. 
% % row vector such as [1 2 3] will be separated by two spaces, that is "1  2  3"
% % One array can contain all of them, but only one value per cell.
% % 2x times faster than Sylvain's codes (8.8s vs. 17.2s):
% % tic;C={'te','tm';5,[1,2];true,{}};C=repmat(C,[10000,1]);cell2csv([datestr(now,'MMSS') '.csv'],C);toc;
% % 
% % Modified and Optimized by Jeremiah Valenzuela, Oct 24, 2019
% % - Improved performance drastically for network drives and a bit for local drives. Faster
% %   than writecell()
% % - Solved escape character issue with the same performance fix
% % - Added ability to use different file permissions, like append
% % - Added ability to specify float precision

%% Checking for optional Variables
if ~exist('delimiter', 'var')
    delimiter = ',';
end

if ~exist('filePermissions', 'var')
    filePermissions = 'w';
end

if ~exist('excelYear', 'var')
    excelYear = 1997;
end

if ~exist('decimal', 'var')
    decimal = '.';
end

%% Setting delimiter for newer excelYears
if excelYear > 2000
    delimiter = ';';
end

% convert cell
cellArray = cellfun(@StringX, cellArray, 'UniformOutput', false);

%% Write file
datei = fopen(fileName, filePermissions);

% Join cell rows with delimiter and convert to string. This drastically
% improves performance of writing the file, especially on network drives,
% and solves any special character issues (like \ and % which are escape
% characters in fprintf.
fprintf(datei, '%s\n', string(join(cellArray, delimiter)));

% Closing file
fclose(datei);

% sub-function
function x = StringX(x)
    % If zero, then empty cell
    if isempty(x)
        x = '';
    % If numeric -> String, e.g. 1, [1 2]
    elseif isnumeric(x) && isrow(x)
        if ~exist('numericPrecision', 'var')
            x = num2str(x);
        else
            if mod(x, 1) % True if decimal
                x = num2str(x, numericPrecision);
            else
                x = num2str(x);
            end
        end
        if decimal ~= '.'
            x = strrep(x, '.', decimal);
        end
    % If logical -> 'true' or 'false'
    elseif islogical(x)
        if x == 1
            x = 'TRUE';
        else
            x = 'FALSE';
        end
    % If matrix array -> a1 a2 a3. e.g. [1 2 3]
    % also catch string or char here
    elseif isrow(x) && ~iscell(x)
        x = num2str(x);
    % everthing else, such as [1;2], {1}
    else
        x = 'NA';
    end

    % If newer version of Excel -> Quotes 4 Strings
    if excelYear > 2000
        x = ['"' x '"'];
    end
end % end sub-function
end % end function