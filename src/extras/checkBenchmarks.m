function checkBenchmarks(directory)
%CHECKBENCHMARKS  check properties of all benchmarks in a directory
%
% Syntax:
%   CHECKBENCHMARKS
%   CHECKBENCHMARKS(directory)
%
% Description:
%   When operating with new benchmark systems for the first time, it is
%   often desired to get an overview of the properties of the systems in 
%   orde to choose the right one for simulations.
%
%   This function loads all .mat files within a directory (so make sure
%   only benchmark files have a .mat extension) and generates a report
%   stating
%   - state-space dimension
%   - number of inputs and outputs
%   - ODE or DAE
%   - stability and strictly dissipativity
%
%
% See also:
%   LOADSSS, SSS\ISSTABLE, SSS\ISSD, CONDEST
%
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  14 Aug 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%   Go to benchmark directory
currDir = cd;
if ~exist('directory','var')
    directory = 'P:\MOR\Matlab\Benchmark\SLICOT';
end
cd(directory);

%   Get all .mat files
benchmarks = dir('*.mat');
fname = '_benchmarkReport.txt';
fid = fopen(fname,'w'); 
% fid = 1; %for debugging
fprintf(fid,'Benchmark report started %s\n',datestr(now));
fprintf(fid,'--------------------------------------------------\n\n');

%   Check all benchmarks and write report

for iBench = 1:length(benchmarks)
    fprintf(fid,'Loading "%s"...\n',benchmarks(iBench).name);
    try
        sys = loadSss(benchmarks(iBench).name);
        info = disp(sys); 
        fprintf(fid,'%s\n',info);

        % check if DAE
        condE = condest(sys.e);
        fprintf(fid,'condE: \t %e\n',condE);

        % check stability and strict dissipativity
        try %stability
            stab = isstable(sys);
            fprintf(fid,'stable: \t%i\n',stab);
        catch err
            warning('isstable failed with message: %s.',err.message);
            fprintf(fid,'isstable failed with message: %s.',err.message);
        end
        try %strict dissipativity
            [sd, mu] = issd(sys);
            fprintf(fid,'issd: \t%i \t\t mu = %e\n',sd,mu);
        catch err
            warning('issd failed with message: %s.',err.message);
            fprintf(fid,'issd failed with message: %s.',err.message);
        end
        fprintf(fid,'\n\n');
    catch err
        fprintf(fid,strrep(getReport(err,'basic','hyperlinks','off'),sprintf('\n'),'-'));
        fprintf(fid,'\n\n');
    end
end
fclose(fid); cd(currDir)

