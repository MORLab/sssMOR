
function x=suggestVarname(x,hObject)

% Funktion schlägt Name für Variable vor, den es noch nicht im workspace
% gibt, angefangen wird mit x, das Vorschlag wird, wenn hObjct übergeben
% wurde, in den String von hObject geschrieben und auch zurück gegeben
% letzte Änderung 11.08.2010
if existInBaseWs(x)==0
    set(hObject,'UserData',0)
    set(hObject,'string',x) 
    return
end
i=2;
while existInBaseWs(sprintf('%s%i',x,i))==1
    i=i+1;
end
x=sprintf('%s%i',x,i);
if isempty(hObject)
    return
end
set(hObject,'UserData',0)
set(hObject,'string',x)
