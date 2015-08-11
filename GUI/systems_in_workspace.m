function x = list_systems_in_workspace()
% finds and lists all dynamical systems that are contained in workspace
    % read all variable names
    s=evalin('base', 'whos');
    % preallocate memory
    x=cell(length(s),1);
    for i=1:length(s)
        if strcmp(s(i).class,'ss') || strcmp(s(i).class,'sss') && ...
            ~strcmp(s(i).name,'load_var_dont_destroy') && ...
            ~strcmp(s(i).name,'GUI_dont_destroy')
            % save name
            x{i}=s(i).name;
        end
    end
    % remove empty (non-system) entries
    x(cellfun(@isempty,x)) = [];