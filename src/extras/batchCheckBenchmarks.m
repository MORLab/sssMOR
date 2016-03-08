%   Define email parameters
setpref('Internet','E_mail','a.castagnotto@tum.de');
setpref('Internet','SMTP_Server','mailout.lrz-muenchen.de');
[~, pcName] = system('hostname'); pcName = strrep(pcName,sprintf('\n'),'');

%   Define directories
% directory = {'P:\MOR\Matlab\Benchmark\SLICOT_test'}; %for debugging

directory = {'P:\MOR\Matlab\Benchmark\SLICOT',...
                'P:\MOR\Matlab\Benchmark\feng',...
                'P:\MOR\Matlab\Benchmark\H-Bridge'};

parentDir = {'P:\MOR\Matlab\Benchmark\IMTEK',...
              'P:\MOR\Matlab\Benchmark\Han',...
              'P:\MOR\Matlab\Benchmark\MORWiki'}; %nested directories
for iParDir = 1:length(parentDir)
    temp = dir(parentDir{iParDir});
    for iDir = 3:length(temp) %the first two are dots
        directory = {directory{:}, fullfile(parentDir{iParDir},temp(iDir).name)};             
    end
end

%   Check Benchmarks
errCount = 0; %start error count
errDir = {};
errMsg = {};
errDetails = '';
try
    tic
    for iDir = 1:length(directory) 
        try
            checkBenchmarks(directory{iDir})
        catch err
            errCount = errCount + 1;
            errDir = {errDir{:}, directory{iDir}};
            errMsg = {errMsg{:}, ...
                strrep(getReport(err,'basic','hyperlinks','off'),sprintf('\n'),'-')};
            errDetails = sprintf('%s\n--\n%s',...
                errDetails,getReport(err,'extended','hyperlinks','off'));
        end
    end
    t = toc;
catch err
    subject = sprintf('Error on %s.',pcName);
    content = sprintf(['Execution stopped during %s with error message:\n%s\n',...
        'The current directory was %s'],...
        mfilename,err.message, directories(iDirs).name);
    sendmail({'a.castagnotto@tum.de', 'a.castagnotto@hotmail.com'},...
        subject,content);
    rethrow(err)
end

%   Generate a formatted string with error messages
errTabFName = 'batchCheckBenchmarks_errorReport.txt';
if errCount == 0
    fid = fopen(errTabFName,'w');
    fprintf(fid,'No errors during the execution!');
    fclose(fid);
else
    errTable = cell2table([errDir', errMsg'],'VariableNames',{'errDir','errMsg'});
    writetable(errTable,errTabFName,'Delimiter',',');
end

%   Conclude execution
subject = sprintf('Execution terminated on %s.',pcName);
content = sprintf(['Execution terminated after %f minutes.\n',...
          '%i directories out of %i stopped prematurely. A summary of the ',...
          'error is attached to this mail.'],...
          t/60, errCount, length(directory));
appendix = errTabFName;
sendmail({'a.castagnotto@tum.de', 'a.castagnotto@hotmail.com'},...
    subject,content,appendix);