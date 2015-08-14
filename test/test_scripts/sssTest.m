function sssTest

%% loadSss

sys = loadSss('build');
sys = loadSss('gyro');

%%  isstable
sys = loadSss('build');
issd(sys), pause
sd = issd(sys), pause
[sd, mu] = issd(sys), pause

sys = loadSss('beam');
issd(sys), pause
[sd, mu] = issd(sys)