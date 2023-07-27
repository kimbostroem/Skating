function saveTable(myTable, fname, fext, fdir)

if nargin < 4
    fdir = '.';
end

if nargin < 3
    fext = '.xlsx';
end

% make file extension lowercase
fext = lower(fext);

% ensure that file extension is cell array
fext = cellstr(fext);

for iFext = 1:length(fext)
    myFext = fext{iFext};
    if ~startsWith(fext, '.')
        myFext = sprintf('.%s', myFext);
    end
    switch myFext
        case '.xlsx'
            % save as Excel table
            outpath = fullfile(fdir, [fname, myFext]);
            writetable(myTable, outpath, 'WriteMode', 'replacefile');
        case '.csv'
            % save as CSV table
            outpath = fullfile(fdir, [fname, myFext]);
            writetable(myTable, outpath, 'WriteMode', 'overwrite', 'Delimiter', ';', 'QuoteStrings', 'all');
        otherwise
            error('Unknown format %s', myFext);
    end
end

end