
function i = exist_in_base_ws(varname)
% returns 1 if x exists in base workspace, otherwise 0
i=evalin('base',sprintf('exist(''%s'')',varname));