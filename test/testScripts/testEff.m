function [] = testEff()
% Evaluates the effectiveness of tbr called with different options
% and input parameters

%choose type 'full' and a high number
loadBenchmarks; 

temp=load('benchmarksSysCell.mat');
sysCell=temp.benchmarksSysCell;
if isempty(sysCell)
    error('No benchmarks loaded.');
end

time1=zeros(length(sysCell),1);
time2=zeros(length(sysCell),1);
time3=zeros(length(sysCell),1);

n=zeros(length(sysCell),1);
n1=zeros(length(sysCell),1);
n2=zeros(length(sysCell),1);
n3=zeros(length(sysCell),1);

% warning off;

sysrCell=cell(length(sysCell),3);

for i=1:length(sysCell)
    sys = sysCell{i};
    if ~sys.isDae && sys.n~=34722 %not peec (dae), not gyro (takes to long)
        disp(['Name: ', sys.Name, 10, 'Order/Inputs/Outputs: ', num2str(sys.n,'%d/'), num2str(sys.m,' %d/'), num2str(sys.p,'%d'), 10,...
            'Descriptor: ', num2str(sys.isDescriptor,'%d')]);
        n(i)=sys.n;
        
        q=150;

        % 2 opts.adi='adi', opts.real=0
        opts.adi='adi';
        opts.real=0;
        opts.sym=0;
        start=tic();
        try
            [sysr2] = tbr(sys, q, opts);
            time2(i)=toc(start);
            n2(i)=sysr2.n;
            sysrCell{i}{2}=sysr2;
            q=sysr2.n;
        catch ME
            time2(i)=Inf;
            n2(i)=0;
            disp(ME.message);
        end
        disp(time2(i));
        disp(n2(i));
        

        % 3 opts.adi='adi', opts.real='real'
        opts.adi='adi';
        opts.real='real';
        opts.sym=0;
        start=tic();
        try
            [sysr3] = tbr(sys, q, opts);
            time3(i)=toc(start);  
            q=sysr3.n;
            n3(i)=sysr3.n;
            sysrCell{i}{3}=sysr3;
        catch ME
            time3(i)=Inf;
            n3(i)=0;
            disp(ME.message);
        end
        disp(time3(i));
        disp(n3(i));
        
        % 1 opts.adi=0
        opts.adi=0;
        opts.real=0;
        opts.sym=0;
        if sys.n~=5177 %not rail_5177 (takes to long)
            start1=tic();
            try
                [sysr1] = tbr(sys, q, opts);
                time1(i)=toc(start1);
                n1(i)=sysr1.n;
                sysrCell{i}{1}=sysr1;
            catch ME
                time1(i)=Inf;
                n1(i)=0;
                disp(ME.message);
            end
        end
        disp(time1(i));
        disp(n1(i));
    end
end

disp([10, 'Result:',10]);

disp([10,'Shortest time: tic-toc']);
disp([num2str(nnz(time1<time2 & time1<time3 & time1>zeros(length(sysCell),1)),'%i') ' 0']);
disp([num2str(nnz(time2<time1 & time2<time3 & time2>zeros(length(sysCell),1)),'%i') ' adi']);
disp([num2str(nnz(time3<time1 & time3<time2 & time3>zeros(length(sysCell),1)),'%i') ' adi + real']);

disp([10,'Reduced order']);
disp([num2str(nnz(n1>=n2 & n1>=n3 & n1>=zeros(length(sysCell),1)),'%i') ' 0']);
disp([num2str(nnz(n2>=n1 & n2>=n3 & n2>=zeros(length(sysCell),1)),'%i') ' adi']);
disp([num2str(nnz(n3>=n1 & n3>=n2 & n3>=zeros(length(sysCell),1)),'%i') ' adi + real']);

delete('benchmarksSysCell.mat');

disp('time:');
disp([time1,time2,time3]);

disp('n');
disp([n1, n2, n3]);
disp(n);

save sysrCellRand;
end