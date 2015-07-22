%% Define version/functions to test

dirname = spritf(...

% go to target directory
cd()

%% Test Bode
%1) bode plot
load(build)
sys = sss(A,B,C,D,E)

bode(sys)

%% Test CURE

load build; %beam
sys = sss(A,B,C);

Opts.CURE.init = 2;
Opts.CURE.test = 1;
Opts.SPARK.test = 1;
Opts.MESPARK.ritz = 1;
sysr = CURE(sys,Opts);


%  Test initializeShifts (within CURE)
for ii = 1:20
    Opts.CURE.nk = ii;
    Opts.CURE.init = 1;
    Opts.zeroThres = 0;
    s0 = initializeShifts(sys,Opts);
    plot(real(s0),imag(s0),'*')
    title(sprintf('Opts.CURE.nk = %i',ii));
    pause
end

