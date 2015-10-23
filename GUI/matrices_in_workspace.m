function x=matrices_in_workspace
% Returns all matrices from workspace

s=evalin('base', 'whos');
% preallocate memory
x=cell(length(s),1); %alle, auch quadtratische
for i=1:length(s)
    if strcmp(s(i).class,'double') && length(s(i).size)==2 && any(s(i).size)
        % save name
        x{i}=s(i).name;
    end
end
% remove empty (non-system) entries
x(cellfun(@isempty,x)) = [];
