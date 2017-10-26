%% Basic test
clc, clear
OptsStruct = struct();
OptsStruct.first     = {'1', '2'};
OptsStruct.second    = {'T', 'F', 'L'};
OptsStruct.third     = {'A', 'B'};

[AllOptsCell,nPermutations] = generateAllOptsNested(OptsStruct);
for i = 1:nPermutations, disp(AllOptsCell{i}), end

%% inner struct test
clc, clear
OptsStruct = struct();
nested.second       = {'T', 'F', 'L'};
nested.third        = {'A', 'B'};

[AllOptsCell,nPermutations] = generateAllOptsNested(nested);
for i = 1:nPermutations, disp(AllOptsCell{i}), end

%% Nested struct test
clc, clear
OptsStruct = struct();
OptsStruct.first               = {'1', '2'};
OptsStruct.nested.second       = {'T', 'F', 'L'};
OptsStruct.nested.third        = {'A', 'B'};

[AllOptsCell,nPermutations] = generateAllOptsNested(OptsStruct);
for i = 1:nPermutations, disp(AllOptsCell{i}), end
