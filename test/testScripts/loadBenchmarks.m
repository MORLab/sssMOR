function [testpath] = loadBenchmarks(Opts)
%% testpath
clear;
p = mfilename('fullpath'); k = strfind(p, '\test'); 
testpath = [p(1:k(end)-1),'\testScripts'];
cd(testpath);

% Default benchmarks (if testScript is not run by test.m)
Def.option = 'light'; %'light','full','heavy'
Def.size =400'; %'light': only benchmarks with sys.n<=Opts.size are tested
                %'heavy': only benchmarks with sys.n>Opts.size are tested
Def.number = 3; %choose maximum number of tested benchmarks

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%% Load benchmarks
%the directory "benchmark" is in sssMOR
p = mfilename('fullpath'); k = strfind(p, 'test\'); 
pathBenchmarks = [p(1:k-1),'benchmarks'];
cd(pathBenchmarks);
badBenchmarks = {'LF10.mat','beam.mat','random.mat',...
    'SpiralInductorPeec.mat'};  

% load files
files = dir('*.mat'); 
benchmarksSysCell=cell(1,Opts.number);
nLoaded=1; %count of loaded benchmarks
disp('Loaded systems:');

warning('off');
for i=1:length(files)
    if nLoaded<Opts.number+1
        switch(Opts.option)
            case 'light'
                sys = loadSss(files(i).name);
                if ~any(strcmp(files(i).name,badBenchmarks)) && size(sys.A,1)<=Opts.size
                    benchmarksSysCell{nLoaded}=sys;
                    nLoaded=nLoaded+1;
                    disp(files(i).name);
                end
            case 'full' 
                benchmarksSysCell{nLoaded} = loadSss(files(i).name);
                nLoaded=nLoaded+1;
                disp(files(i).name);
            case 'heavy'
                sys = loadSss(files(i).name);
                if any(strcmp(files(i).name,badBenchmarks)) || size(sys.A,1)>Opts.size
                    benchmarksSysCell{nLoaded}=sys;
                    nLoaded=nLoaded+1;
                    disp(files(i).name);
                end
            otherwise
                error('Benchmark option is wrong.');
        end
    end
end
benchmarksSysCell(nLoaded:end)=[];
warning('on');

% change path back and save loaded systems
cd(testpath);
save('benchmarksSysCell.mat');
end