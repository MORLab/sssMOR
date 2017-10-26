function [AllOptsCellArray,nPermutations] = generateAllOptsNested(AllOptsStruct)
    %  Extract fields and values from AllOpts
    [fieldNames, values_CellArray] = getFieldsAndValues(AllOptsStruct);

    % Generate all combinations using "allcomb"
    [Permutations, nPermutations] = getFieldsPermutations(values_CellArray);
    
    % Construct AllOptsCellArray of Opts structures
    AllOptsCellArray = getCellArray(Permutations, nPermutations, fieldNames);
end

function [fieldNames, values_CellArray] = getFieldsAndValues(AllOptsStruct)
    fieldNames = fieldnames(AllOptsStruct).';
    values_CellArray = cell(1,length(fieldNames));

    for iField = 1:length(fieldNames)
        if isstruct(AllOptsStruct.(fieldNames{iField}))
            % Recurcively generate cell array for the nested struct
            values_CellArray{iField} = generateAllOptsNested(AllOptsStruct.(fieldNames{iField}));
        else
            values_CellArray{iField} = AllOptsStruct.(fieldNames{iField});
        end
    end
end

function [Permutations, nPermutations] = getFieldsPermutations(values_CellArray)
    % MATLAB central file exchange
    % https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin-
    % downloaded 23 Nov 2016
    Permutations = allcomb(values_CellArray{:});
    nPermutations = size(Permutations,1);
end

function AllOptsCellArray = getCellArray(Permutations, nPermutations, fieldNames)
    AllOptsCellArray = cell(1,nPermutations);
    for iPermutations = 1:nPermutations
        AllOptsCellArray{iPermutations} = cell2struct(Permutations(iPermutations,:), fieldNames, 2);
    end
end