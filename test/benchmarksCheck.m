function benchmarksCheck(varargin)
%% List of standard benchmarks that should always be in the Benchmarks folder
defaultBenchmarksList = {...
	'CDplayer.mat';...
    'LF10.mat';...
    'SpiralInductorPeec.mat';...
    'beam.mat';...
    'building.mat';...
    'eady.mat';...
    'fom.mat';...
    'gyro.mat';...
    'heat-cont.mat';...
    'iss.mat';...
    'peec.mat';...
    'rail_1357.mat';...
    'rail_5177.mat';...
    'random.mat';...
	};

%% Path of local benchmarks folder
[testFolderPath,~,~] = fileparts(mfilename('fullpath'));
sssmorFolderPath = fileparts(testFolderPath);
benchmarksFolderPath = [sssmorFolderPath filesep 'benchmarks' filesep];

%% Path to the NAS benchmarks folder
nasFolderPath = [filesep filesep 'nas.ads.mwn.de' filesep 'tumw' filesep 'rtk' filesep 'MORLab' filesep 'Toolbox' filesep '05_benchmark' filesep'];

%% Parsing
minInputs = 0;
maxInputs = 2;
narginchk(minInputs,maxInputs)

p = inputParser;
checkPath = @(x) any(exist(x,'dir'));

addOptional(p,'destinationPath',benchmarksFolderPath,checkPath);
addOptional(p,'sourcePath',nasFolderPath,checkPath);
addOptional(p,'benchmarksList',defaultBenchmarksList,@iscellstr);

p.FunctionName = 'benchmarks_check';

parse(p,varargin{:});

%% Shorten variable names
sourcePath = p.Results.sourcePath;
destinationPath = p.Results.destinationPath;
benchmarksList = p.Results.benchmarksList;


%% Determine which mat-files are already in the local folder
filesInFolder = what(destinationPath); % struct with files of all types
matFilesInFolder = filesInFolder.mat; % cell array with names of mat-files

%% Check if all the standard benchmarks are in the local folder
nBenchmarks = length(benchmarksList);

if isempty(matFilesInFolder) % if there are no mat-files in the local folder, then copy all benchmarks
	for iBenchmarks = 1:nBenchmarks
		fprintf(2,'Benchmark %s not found. Copying it...\n...from %s\n...to %s\n\n',...
			benchmarksList{iBenchmarks}, sourcePath, destinationPath);
		copyfile([sourcePath benchmarksList{iBenchmarks}] , destinationPath);
	end
else
	for iBenchmarks = 1:nBenchmarks
		if ~any(strcmp(p.Results.benchmarksList{iBenchmarks}, matFilesInFolder))
			fprintf(2,'Benchmark %s not found. Copying it...\n...from %s\n...to %s\n\n',...
				benchmarksList{iBenchmarks}, sourcePath, destinationPath);
			copyfile([sourcePath benchmarksList{iBenchmarks}] , destinationPath);
		end
	end
end






