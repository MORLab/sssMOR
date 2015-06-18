%% Define version/functions to test

version = 1;

dirname = spritf(...

% go to target directory
cd()

%% Test Bode
%1) bode plot
load(build)
sys = sss(A,B,C,D,E)

bode(sys)

%%