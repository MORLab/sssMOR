function updateIndexFiles( helpSource_, productName_, footer_, doc_dir )
% read the existing helptoc file, check that the processed files match the
% existing in the helptoc and modify it accordingly
%
% 
% The current helptoc.xml backed up in helptoc_backup.xml and is used to
% generate the categories. Hence, to change the category of an entry and/or
% the ordering of a category, only the helptoc.xml has to be changed
% accordingly. When the help is rebuilt the changes will be preserved.
% 
% 
%--------------------------------------------------------------------------
% Copyright (c) ITM University of Stuttgart, <a href="matlab:
%  web('http://www.itm.uni-stuttgart.de')">www.itm.uni-stuttgart.de</a>
%--------------------------------------------------------------------------

%% ChangeLog
%{
---------------------------------------------------------------------------
 #      Date      Name      Description/Comment
===========================================================================
[001]   12/2012   N. Walz   First Version.
---------------------------------------------------------------------------
%}
% END ChangeLog
%% Function Body

helptoc_dir = [doc_dir filesep 'html' filesep 'helptoc.xml'];

% check whether current helptoc file exists
if exist(helptoc_dir, 'file') ~= 2
    % RODRIGO: modified from 'helptoc.xml' to 'doc_dir/html/helptoc.xml', because 
    % builddocsearchdb only works if helptoc.xml is in the folder where the
    % html files are

    error('UPDATEINDEXFILES:FileNonexistend', ...
        ['Help TOC can not be updated since no ''helptoc.xml'' file has been found ' ...
        'in the current directory.\nPlease generate this file first.'])
end

% read helptoc document
% RODRIGO: xDoc = document node. xDoc.getNodeType == 9
xDoc = xmlread( helptoc_dir );
% RODRIGO: modified from 'helptoc.xml' to 'doc_dir/html/helptoc.xml' (read previous
% comment)

% -------------------------------------------------------------------------
% Note: Lists returned by DOM methods use zero-based indexing.
% -------------------------------------------------------------------------

year_ = datestr(now, 'yyyy');
comment_ = sprintf( ...
  [ '\n%s Toolbox\n' ...
  'Copyright (c) %s' ...
  ' Institute of Automatic Control,\n' ...
  '  Technische Universitaet Muenchen\n' ...
  'See LICENSE.txt for license information.\n' ...
  'web: https://www.rt.mw.tum.de/\n' ...
  '\n(This file has been generated automatically)\n' ], ...
  upper(productName_), year_);

% modify the comment
cmtElt_ = xDoc.getFirstChild;
% RODRIGO: cmtElt_ = "header" comment, therefore it is the 
% first child of the document node.
% RODRIGO: cmtElt_ = comment node. cmtElt_.getNodeType == 8

if ~isempty( cmtElt_ ) && strcmp( cmtElt_.getNodeName, '#comment' )
    cmtElt_.setNodeValue( comment_ );
end

toc_ = xDoc.getDocumentElement;
% RODRIGO: toc_ = element node. toc_.getNodeType == 1
% RODRIGO: getDocumentElement method gets the root element node
% of the xDoc document node.



% get the main <tocitem> element, this contains the main product page

% RODRIGO:  remove all empty text nodes that are child of
% the root element node <toc>
l_ = 0;
while l_ < toc_.getLength
    tocItem_ = toc_.item(l_);
    if strcmp(char(tocItem_.getNodeName), '#text') && ...
            isempty(strtrim(char(tocItem_.getNodeValue)))
        toc_.removeChild( tocItem_ );
    else
        l_ = l_ + 1;
    end
end


tocItemList = toc_.getElementsByTagName('tocitem').item(0);
% RODRIGO: getElementsByTagName returns a NodeList of all the 
% element nodes in document order with a given tag name. Each
% element node in the NodeList also includes its child nodes.
% toc_.getElementsByTagName('tocitem') is a NodeList of all 
% descendant elements of <toc> with the <tocitem> tag name, 
% in document order. item(0) of that NodeList is the first 
% (and main) <tocitem> element node. Therefore, tocItemList is 
% the main <tocitem> element node (all other <tocitem> elements 
% are child nodes of it), together with all of its child nodes.



% RODRIGO: The following if-condition checks if the productName_ string
% passed to updateIndexFiles is also in the text element of the main 
% <tocitem> element of the helptoc.xml file. If it isn't it shows a warning.
if isempty( strfind( ...
        lower(char(tocItemList.getFirstChild.getNodeValue)), ...
        lower(productName_) ) )
    warning(['helptoc file does not seem to belong to ' productName_])
end


% RODRIGO:  remove all empty text nodes that are child nodes of the
% main <tocitem> element.
l_ = 0;
while l_ < tocItemList.getLength
    tocItem_ = tocItemList.item(l_);
    if strcmp(char(tocItem_.getNodeName), '#text') && ...
            isempty(strtrim(char(tocItem_.getNodeValue)))
        tocItemList.removeChild( tocItem_ );
    else
        l_ = l_ + 1;
    end
end
% RODRIGO: the main <tocitem> element ('tocItemList') has now non-empty text
% nodes, comment nodes, <tocitem> element nodes and maybe other element nodes




%{  OLD ORIGINAL COMMENTS

% remove trailing whitespace and newline characters
% tocItemList.getFirstChild.setNodeValue( sprintf('%s\n', strtrim(char( ...
%     tocItemList.getFirstChild.getNodeValue))) )

% first level of toc items, 
% rootChilds_ = rootTocItem.getChildNodes;
% tocItemList = rootChilds_.getElementsByTagName( 'tocitem' );
% tocItemList = rootTocItem.getElementsByTagName( 'tocitem' );


% for k_ = 0:tocItemList.getLength-1

% start with second item, since first item has already been processed
k_ = 0;%1;

%}






k_ = 0;
while k_ < tocItemList.getLength
    
    % get a first level entry of the tocItemList nodelist
    thisListItem = tocItemList.item(k_);
    
    % check whether item is a <tocitem> entry
    if ~thisListItem.hasChildNodes 
    % RODRIGO: filter out nodes without child nodes
    
        % delete leading and trailing whitespace and newline characters
        nodeVal_ = strtrim(char( thisListItem.getNodeValue ));
        thisListItem.setNodeValue( sprintf('%s\n', nodeVal_) );
        if isempty( char(thisListItem.getNodeValue) )
            [~] = tocItemList.removeChild( thisListItem );
            % RODRIGO: remove empty nodes (text,tocitems,etc.)
        else
            % increment index
            k_ = k_ + 1;
        end
        continue; % RODRIGO continue to next iteration of main while-loop
    elseif ~strcmp('tocitem', char(thisListItem.getTagName) ) 
    % RODRIGO: filter out nodes with child nodes which 
    % arent <tocitem> element nodes

        % increment index
        k_ = k_ + 1;
        continue; % RODRIGO continue to next iteration of main while-loop
    end
    
    
    % RODRIGO: only <tocitem> element nodes with child nodes make it this
    % far in the while-loop.
    % RODRIGO: remove empty child nodes of the <tocitem> element node
    l_ = 0;
    while l_ < thisListItem.getLength
        tocItem_ = thisListItem.item(l_);
        if strcmp(char(tocItem_.getNodeName), '#text') && ...
                isempty(strtrim(char(tocItem_.getNodeValue)))
            thisListItem.removeChild( tocItem_ );
        else
            l_ = l_ + 1;
        end
    end

    % get the current branch
    % RODRIGO: This is a switch case statement with 4 cases:
    % 'getting started'
    % 'functions reference'
    % 'examples'
    % and the standard otherwise case
    switch lower( strtrim( char( thisListItem.getFirstChild.getNodeValue ) ) )
        
        case 'getting started'
            
            chldItems_ = thisListItem.getElementsByTagName( 'tocitem' );
            
            list_ = {};
            for l_ = 0:chldItems_.getLength-1
                
                % get the target of the current entry
                match = regexp(char( chldItems_.item(l_).getAttribute('target') ), ...
                    '(?=[./]*)([^./]*)(?=\.html)', 'match');
                    % RODRIGO: \w+(?=\.html)
                    % RODRIGO: TEMP = char( chldItems_.item(l_).getAttribute('target') )
                    % [~,match,~] = fileparts(TEMP)
                % RODRIGO: the regexp matches the name of the target html-files
                % under the <tocitem> element named "Getting Started"
                
                
                list_ = [list_; match]; %#ok
            end
            
            % check if help file list matches the tocitem elements
            check_idx = false;
            
            % RODRIGO: helpSource_ is a struct containing information about
            % the marked-up m-files that are going to be published. 
            % Information like category, path and sometime a short 
            % description of the marked-up m-files.
            for iItem_ = 1:length(helpSource_.getting_started.item)
                
                curItem_ = helpSource_.getting_started.item(iItem_);
                [~, filename, ~] = fileparts(curItem_.file);
                
                idx = strcmp(filename, list_);
                check_idx = check_idx | idx;
                
                % create non-existent tocitem elements
                if isempty( find(idx, 1) )
                    % RODRIGO: create new <tocitem> element node
                    newNode_ = xDoc.createElement('tocitem');
                     
                    % RODRIGO: This was the old configuration, 
                    % now that helptoc.xml is in the html folder
                    % the path is different
                    % newNode_.setAttribute( 'target', fullfile('.', ...
                    %   'html', [filename '.html']) );
                    
                    % RODRIGO: New path configuration
                    newNode_.setAttribute( 'target', fullfile([filename '.html']) );

                    if isempty( curItem_.name )
                        txt_ = filename;
                    else
                        txt_ = curItem_.name;
                    end

                    % RODRIGO: append text node to the new <tocitem> element node
                    newNode_.appendChild( xDoc.createTextNode( txt_ ) );
                    % RODRIGO: append new <tocitem> element node to 
                    % "Getting Started" <tocitem> element node
                    
                    thisListItem.appendChild( newNode_ );
                end
            end
            
            % delete elements for which no help file exists
            if ~all(check_idx)
                idx = find(~check_idx);
                for i_ = 1:length(idx)
                    [~] = thisListItem.removeChild( chldItems_.item(idx(i_)) );
                    idx = idx - 1;
                end
            end
            
        case 'functions reference'
            
            % category list
            catList_ = {};
            % content struct
            cont = struct('name', {}, 'file', {}, 'cont', {});
            
            funIndexHtmlFile_ = char(thisListItem.getAttribute('target'));
            if isempty(funIndexHtmlFile_)
                
                funIndexFile_ = fullfile(doc_dir, filesep, 'source', 'functions_index.m');
                
                % RODRIGO: This was the old configuration, 
                % now that helptoc.xml is in the html folder
                % the path is different
                % funIndexHtmlFile_ = fullfile( '.', 'html', 'functions_index.html');
                
                % RODRIGO: New path configuration
                funIndexHtmlFile_ = fullfile('functions_index.html');
                
                thisListItem.setAttribute('target',  funIndexHtmlFile_);
            else
                [~, filename, ~] = fileparts(funIndexHtmlFile_);
                funIndexFile_ = fullfile(doc_dir, 'source', [filename '.m']);
            end
            
            % loop over child elements
            l_ = 1;
            while l_ < thisListItem.getLength
                
                % get the function categories
                % RODRIGO: first sublevel of <tocitem> elements under the
                % "functions reference" <tocitem> are the function categories
                catItem = thisListItem.item(l_);
                
                % check whether item is a tocitem entry
                if ~catItem.hasChildNodes
                    % delete leading and trailing whitespace and newline characters
                    catItem.setNodeValue( strtrim(char( ...
                        catItem.getNodeValue )));
                    if isempty( char(catItem.getNodeValue) )
                        [~] = thisListItem.removeChild( catItem );
                    else
                        % increment index
                        l_ = l_ + 1;
                    end
                    continue;
                elseif ~strcmp('tocitem', char(catItem.getTagName) )
                    % increment index
                    l_ = l_ + 1;
                    continue;
                end
                
                % get category name
                catName_ = strtrim( char( catItem.getFirstChild.getNodeValue ) );
                
                % remove white space in node value
                catItem.getFirstChild.setNodeValue(sprintf('%s\n', catName_));
                
                % store name and category item
                catList_ = [catList_ ; {catName_, catItem } ]; %#ok
                
                % update content struct for index file generation
                cont(l_).name = catName_;
                
                % ensure that the category points to the functions index
                % html page
                catItem.setAttribute('target', funIndexHtmlFile_);
                
                
                l_ = l_ + 1;
            end
            
            % include the 'not yet categorized' category if not present
            if ~any( strcmp('Not Categorized', catList_(:,1)) )
                newNode_ = xDoc.createElement('tocitem');
                newNode_.setAttribute('target', funIndexHtmlFile_);
                newNode_.appendChild( xDoc.createTextNode('Not Categorized') );
                thisListItem.appendChild( newNode_ );
                % add to catList_
                catList_ = [catList_; {'Not Categorized', newNode_}]; %#ok
                cont(end+1).name = 'Not Categorized';
            end
           
            % check if help file list matches the tocitem elements
            check_idx = 1:length(helpSource_.functions.item);
            
            for iC_ = 1:size(catList_,1)
                
                funItemList_ = catList_{iC_,2};
                
                iItem_ = 1;
                while iItem_ < funItemList_.getLength
                    
                    curFunItem_ = funItemList_.item( iItem_ );
                    
                    % check whether item is a tocitem entry
                    if ~curFunItem_.hasChildNodes
                        % delete leading and trailing whitespace and newline characters
                        curFunItem_.setNodeValue( strtrim(char( ...
                            curFunItem_.getNodeValue )));
                        % if is empty text node, remove it
                        if isempty( char(curFunItem_.getNodeValue) )
                            [~] = funItemList_.removeChild( curFunItem_ );
                        else
                            % increment index
                            iItem_ = iItem_ + 1;
                        end
                        continue;
                    elseif ~strcmp('tocitem', char(curFunItem_.getTagName) )
                        % increment index
                        iItem_ = iItem_ + 1;
                        continue;
                    end
                    
                    
                    curFunName_ = strtrim( char(curFunItem_.getFirstChild.getNodeValue) );
                    
                    idx = strcmp( {helpSource_.functions.item(check_idx).name}, ...
                         curFunName_ );
                
                    
                    % delete elements for which no help file exists
                    if ~any(idx)
                        [~] = catList_{iC_,2}.removeChild( curFunItem_ );
                    else
                        [~, filename, ~] = fileparts( ...
                            helpSource_.functions.item(check_idx(idx)).file );

                        % RODRIGO: This was the old configuration, 
                        % now that helptoc.xml is in the html folder
                        % the path is different
                        % curFunItem_.setAttribute('target', fullfile('.', ...
                        %    'html', [filename '.html']));
                        
                        % RODRIGO: New path configuration
                        curFunItem_.setAttribute( 'target', fullfile([filename '.html']) );

                        curFunItem_.setAttribute('image', 'HelpIcon.FUNCTION');
                        
                        cont(iC_).cont(end+1).name = curFunName_;
                        
                        % RODRIGO: This path doesn't need to be changed
                        % because "functions_index.hmtl" is in the same
                        % folder as all the other html-files
                        cont(iC_).cont(end).file = fullfile('.', ...
                            [filename '.html']);

                        cont(iC_).cont(end).description = ...
                            helpSource_.functions.item(check_idx(idx)).description;
                        
                        iItem_ = iItem_ + 1;
                    end
                    
                    % remove idx from check list, if found any match
                    check_idx(idx) = [];
                    
                end
                
            end
            
            
            % create elements for non-existent tocitem elements
            if ~isempty( check_idx )
                uncatIdx_ = strcmp('Not Categorized', catList_(:,1));
                for iC_ = check_idx
                    curItem_ = helpSource_.functions.item(iC_);
                    [~, filename, ~] = fileparts(curItem_.file);
                    newNode_ = xDoc.createElement('tocitem');

                    % RODRIGO: This was the old configuration, 
                    % now that helptoc.xml is in the html folder
                    % the path is different
                    % newNode_.setAttribute( 'target', fullfile( ...
                    %    '.', 'html', [filename '.html']) );

                    % RODRIGO: New path configuration
                    newNode_.setAttribute( 'target', fullfile([filename '.html']) );

                    newNode_.setAttribute('image', 'HelpIcon.FUNCTION');
                    if isempty( curItem_.name )
                        txt_ = filename;
                    else
                        txt_ = curItem_.name;
                    end
                    newNode_.appendChild( xDoc.createTextNode( txt_ ) );
                    catList_{uncatIdx_,2}.appendChild( newNode_ );
                    
                    cont(uncatIdx_).cont(end+1).name = txt_; %curFunName_;

                    % RODRIGO: This path doesn't need to be changed
                    % because "functions_index.hmtl" is in the same
                    % folder as all the other html-files
                    cont(uncatIdx_).cont(end).file = fullfile('.', ...
                        [filename '.html']);

                    cont(uncatIdx_).cont(end).description = curItem_.description;
                end
            end
            
            % remove category entry if it is empty
            for iC_ = 1:size(catList_,1)
                if catList_{iC_,2}.getElementsByTagName('tocitem').getLength == 0
                    [~] = thisListItem.removeChild( catList_{iC_,2} );
                    cont(iC_) = [];
                end
            end
            
            
            % write the functions index file
            writeIndexFile( funIndexFile_, 'Functions Reference', ...
                ['UI-functions of the ' productName_ ' toolbox sorted by categories'], ...
                cont);
            
        case 'examples'
            
            chldItems_ = thisListItem.getElementsByTagName( 'tocitem' );
            
            list_ = {};
            for l_ = 0:chldItems_.getLength-1
                
                % get the target of the current entry
                match = regexp(char( chldItems_.item(l_).getAttribute('target') ), ...
                    '(?=[./]*)([^./]*)(?=\.html)', 'match');
                    % RODRIGO: \w+(?=\.html)
                    % RODRIGO: TEMP = char( chldItems_.item(l_).getAttribute('target') )
                    % [~,match,~] = fileparts(TEMP)
                
                list_ = [list_; match]; %#ok
            end
            
            % check if help file list matches the tocitem elements
            check_idx = false;
            for iItem_ = 1:length(helpSource_.examples.item)
                
                curItem_ = helpSource_.examples.item(iItem_);
                [~, filename, ~] = fileparts(curItem_.file);
                
                idx = strcmp(filename, list_);
                check_idx = check_idx | idx;
                
                % create non-existent tocitem elements
                if isempty( find(idx, 1) )
                    newNode_ = xDoc.createElement('tocitem');

                    % RODRIGO: This was the old configuration, 
                    % now that helptoc.xml is in the html folder
                    % the path is different
                    % newNode_.setAttribute( 'target', fullfile('.', ...
                    %   'html', [filename '.html']) );
                    
                    % RODRIGO: New path configuration
                    newNode_.setAttribute( 'target', fullfile([filename '.html']) );

                    if isempty( curItem_.name )
                        txt_ = filename;
                    else
                        txt_ = curItem_.name;
                    end
                    newNode_.appendChild( xDoc.createTextNode( txt_ ) );
                    thisListItem.appendChild( newNode_ );
                end
            end
            
            % delete elements for which no help file exists
            if ~all(check_idx) && ~isempty(helpSource_.examples.item)
                idx = find(~check_idx);
                for i_ = 1:length(idx)
                    [~] = thisListItem.removeChild( chldItems_.item(idx(i_)-1) );
                    idx = idx - 1;
                end
            end
            
        otherwise
            % leave as is
    end
    
    % increment index
    k_ = k_ + 1;
end


xmlwrite(helptoc_dir, xDoc);

% str_ = xmlwrite(xDoc);
% disp(str_)

end


%% subfunctions
% RODRIGO: This subfunction writes the "functions_index.m" file
function writeIndexFile(filename, head, intro, cont)


str_ = sprintf('%%%% %s\n', head);
str_ = sprintf('%s%% %s\n', str_, intro);

for iC_ = 1:length(cont)
    str_ = sprintf('%s\n%%%% %s\n%%\n', str_, strtrim(cont(iC_).name));
    for iE_ = 1:length(cont(iC_).cont)
        str_ = sprintf('%s%% * <%s |%s|>  --  %s\n', str_, ...
            cont(iC_).cont(iE_).file, cont(iC_).cont(iE_).name, ...
            cont(iC_).cont(iE_).description);
    end
end

% RODRIGO: take the footer defined by publishHelp
strFooter_ = evalin('caller', 'footer_');
str_ = sprintf('%s%%\n%s', str_, strFooter_);
% str_ = sprintf('%s%%\n%s', str_, getFooter);

fid = fopen(filename, 'w');
fprintf(fid, '%s', str_);
fclose(fid);

end


function writeExamplesIndexFile(filename, head, intro, cont)


str_ = sprintf('%%%% %s\n', head);
str_ = sprintf('%s%% %s\n', str_, intro);

for iC_ = 1:length(cont)
    str_ = sprintf('%s\n%%%% %s\n%%\n', str_, strtrim(cont(iC_).name));
    for iE_ = 1:length(cont(iC_).cont)
        str_ = sprintf('%s%% * <%s |%s|>  --  %s\n', str_, ...
            cont(iC_).cont(iE_).file, cont(iC_).cont(iE_).name, ...
            cont(iC_).cont(iE_).description);
    end
end

% RODRIGO: take the footer defined by publishHelp
strFooter_ = evalin('caller', 'footer_');
str_ = sprintf('%s%%\n%s', str_, strFooter_);
% str_ = sprintf('%s%%\n%s', str_, getFooter);

fid = fopen(filename, 'w');
fprintf(fid, '%s', str_);
fclose(fid);

end


function writeInfoIndexFile(filename, head, intro, cont)


str_ = sprintf('%%%% %s\n', head);
str_ = sprintf('%s%% %s\n', str_, intro);

for iC_ = 1:length(cont)
    str_ = sprintf('%s\n%%%% %s\n%%\n', str_, strtrim(cont(iC_).name));
    for iE_ = 1:length(cont(iC_).cont)
        str_ = sprintf('%s%% * <%s |%s|>  --  %s\n', str_, ...
            cont(iC_).cont(iE_).file, cont(iC_).cont(iE_).name, ...
            cont(iC_).cont(iE_).description);
    end
end

% RODRIGO: take the footer defined by publishHelp
strFooter_ = evalin('caller', 'footer_');
str_ = sprintf('%s%%\n%s', str_, strFooter_);
% str_ = sprintf('%s%%\n%s', str_, getFooter);

fid = fopen(filename, 'w');
fprintf(fid, '%s', str_);
fclose(fid);

end

% function str_ = getFooter
% 
% year_ = datestr(now, 'yyyy');
% str_ = sprintf( ...
%     ['%%%%\n', ...
%     '%% <html>\n', ...
%     '%%   <hr>\n', ...
%     '%%   <p class="copy">&copy; %s ITM University of Stuttgart\n', ...
%     '%%        <tt class="minicdot">&#149;</tt>\n', ...
%     '%%        <a href="http://www.itm.uni-stuttgart.de">Website</a>\n', ...
%     '%%        <tt class="minicdot">&#149;</tt>\n', ...
%     '%%        <a href="file:./LICENSE.html">Terms of Use</a>\n', ...
%     '%%   </p>\n', ...
%     '%% </html>\n'], year_);
% 
% end