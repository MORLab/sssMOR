function [] = compareOrthogonalization()
% Compares the effectiveness of reorthogonalization with dkgs or gs by 
% evaluating the orthogonality of V and the required time.

OrthPoints1=0;
OrthPoints2=0;
TimePoints1=0;
TimePoints2=0;

loadBenchmarks; %choose type 'full'

temp=load('benchmarksSysCell.mat');
sysCell=temp.benchmarksSysCell;
if isempty(sysCell)
    error('No benchmarks loaded.');
end

for i=1:length(sysCell)
    sys = sysCell{i};
    sys = sys(1,1); E = sys.E; A = sys.A; B = sys.B;
    opts.dkgs=1;
    opts.reorth=0;
    start1=tic();
    [V1] = arnoldi(E,A,B,[Inf, 0, 100, 4+13i, 4-13i], @(x,y) (x'*y), opts);
    time1=toc(start1);
    opts.dkgs=0;
    opts.reorth='gs';
    start2=tic();
    [V2] = arnoldi(E,A,B,[Inf, 0, 100, 4+13i, 4-13i], @(x,y) (x'*y), opts);
    time2=toc(start2);
    normV1=norm(V1'*V1-speye(size(V1,2)),'fro');
    normV2=norm(V2'*V2-speye(size(V1,2)),'fro');
    disp([10, sysCell{i}.Name]);
    if normV1==normV2
        OrthPoints1=OrthPoints1+1;
        OrthPoints2=OrthPoints2+1;
        disp('Orth: dkgs and Reorth');
    elseif normV1>normV2
        OrthPoints2=OrthPoints2+1;
        disp('Orth: Reorth');
    else
        OrthPoints1=OrthPoints1+1;
        disp('Orth: dkgs');
    end
    if time1==time2
        TimePoints1=TimePoints1+1;
        TimePoints2=TimePoints2+1;
        disp('Time: dkgs und Reorth');
    elseif time1>time2
        TimePoints2=TimePoints2+1;
        disp('Time: Reorth');
    else
        TimePoints1=TimePoints1+1;
        disp('Time: dkgs');
    end
end

disp([10, 'Result:']);
disp(['Best orthgonality: ' num2str(OrthPoints1,'%i') ' dkgs - ' num2str(OrthPoints2,'%i') ' Reorth']);
disp(['Fastest time: ' num2str(TimePoints1,'%i') ' dkgs - ' num2str(TimePoints2,'%i') ' Reorth']);


delete('benchmarksSysCell.mat');
end