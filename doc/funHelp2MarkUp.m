function funRef_ = funHelp2MarkUp(varargin)
%funHelp2MarkUp  extract function help comments and translate to mark-up
%                for publish function
% Syntax:
%   funHelp2MarkUp('prop', val, ...)
%
% Input Arguments:
% 	-'Files': Cell array of strings; 
% 	-'OutputDir':, String {'./'};
% 	-'Footer': String {''}; If a non-empty string is given it will be
%           appended to every generated help m-file.
% 	-'Format': String ({'html'}, 'latex'); Output format.
% 	-'NSpace': Numeric scalar {4}; Number of spaces to be used for tab
%           replacement.
% 	-'WriteUnmatched: Boolean {true}; Toggles whether unmatched content
%           is written unmodifed, i.e. without markup, or is dismissed from
%           the output.
% 	-'RepFunNames': Boolean {false}; If set to true, not only the name of
%           the current function is marked up to be displayed in typewriter
%           font, but all functions that are given by 'Files' are used.
%           I.e. function names are marked up across all given help files.
%   -'FormatErrors': Boolean {false}; If set to true, non-conforming format
%           of a help comment issues an error. Otherwise, a warning is
%           displayed and the file is skipped.
%   -'OutSuffix': Character string {'help'}; Suffix of the output
%           file containing the marked-up help text, i.e. markup of 
%           'fun_file.ext' is written to 'fun_filehelp.ext' using the
%           default value.
%   -'SkipPattern': Cell array of character strings {{} empty cell}; 
%           Pattern given as regular expression that when encountered in
%           the parsed m-file causes the line to be skipped from
%           processing. It is possible to provide several skip patterns.
%   -'StopPattern': Cell array of character strings {{'^%-{5,}$'}}; 
%           Pattern given as regular expression that when encountered in
%           the parsed m-file terminates the markup process. This is useful
%           in order to skip footer or signature lines in the help
%           comments. It is possible to provide several skip patterns.
%   -'RegExp': Struct {[]}; May be used to modify the regular expression for
%           the markup translation.
%   -'TableBlankCorrection': Numeric {1}; The amount of blanks that is
%           added to the blank counted by the regexp pattern for tables. As
%           tables and table cells are characterized by identation levels
%           this correction can be used to specify if a following line is
%           part of the table / cell or not. This might be necessary for
%           some table patterns.
% 
% Output Arguments:
%   -funRef_: filenames of the created and marked-up help m-files
%           including the output directory.
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
[001]   02/2012   N. Walz   First Version.
[002]   06/2013   N. Walz   Removed dependecy on FAMOUS, so that the
                            function is usable with other programs as well.
---------------------------------------------------------------------------
%}
% END ChangeLog
%% Function Body

% -------------------------------------------------------------------------
% regular expression for matching the

% heading
regExpPat.head = '^% (?<head>[^\s].*):\s*'; % heading
% table 
regExpPat.table = ['^%(?<nblank>\s{2,})' ... table
                '-(?<prop>[^:]+):' ... property
                '\s*(?<text>[^\s]?.*)$']; % cell content
% note box
regExpPat.info = ['^%(?<nblank>\s*)' ...
                '//\s*(?<info_title>[^\s][^:]*):?' ... title line
                '\s*(?<info_text>[^\s:][^:]*)$']; % content of title line
% source code
regExpPat.source = ['^%(?<nblank>\s*)' ...
                '>(?<source>.*)$']; % source code
% plain text
regExpPat.text = ['^%(?<nblank>\s{2,})' ...
                '(?<text>[^\s].*)$']; % plain text
% -------------------------------------------------------------------------
            

% parse input arguments
op = inputParser;

op.addParamValue('Files', {}, @iscellstr);
op.addParamValue('OutputDir', './', @ischar);
op.addParamValue('Footer', '', @ischar);
op.addParamValue('Format', 'html', @(x) any(strcmp({'html', 'latex'}, x)) );
op.addParamValue('NSpace', 4, @isscalar);
op.addParamValue('WriteUnmatched', true, @islogical);
op.addParamValue('RepFunNames', false, @islogical);
op.addParamValue('FormatErrors', false, @islogical);
op.addParamValue('OutSuffix', 'help', @ischar);
op.addParamValue('SkipPattern', {}, @iscellstr);
op.addParamValue('StopPattern', {'^%-{5,}\s*$'}, @iscellstr);
op.addParamValue('RegExp', [], @isstruct)
op.addParamValue('TableBlankCorrection', 1, @isnumeric);


op.parse( varargin{:} );

opts = op.Results;

% modify regular expression
if ~isempty(opts.RegExp)
    fields_ = fieldnames(opts.RegExp);
    for iF_ = 1:length(fields_)
        regExpPat.(fields_{iF_}) = opts.RegExp.(fields_{iF_});
    end
end


% the files to be parsed
help_files = opts.Files;

% check for existence of output directory
if exist(opts.OutputDir, 'dir') ~= 7
    mkdir( opts.OutputDir );
end

% ensure absolute paths
cur_dir = pwd;  % save current path
cd( cur_dir );
cd( opts.OutputDir );
OutputDir = pwd;
cd( cur_dir );  % return to origin


% repFunNames_ contains all function names that are currently processed
% if option is true, these names are marked up across all processed files
repFunNames_ = {};
for iF_ = 1:length(help_files)
    [~, repFunNames_{iF_}, ~] = fileparts(help_files{iF_}); %#ok
end


fprintf('Parsing m-files to extract the function help.\n');
fprintf('\tOutput directory:\n\t\t%s\n', OutputDir);

iRef_ = 0;
for f_ = 1:length(help_files)
    
    fprintf('Current File:\n\t''%s''\n', help_files{f_});
    
    [help_files{f_}, desc_] = parseMFile( help_files{f_}, OutputDir);
    
    fprintf('... done.\n')
    
    if ~isempty(help_files{f_})
        iRef_ = iRef_ + 1;
        funRef_(iRef_).file = [repFunNames_{f_} opts.OutSuffix '.m']; %#ok
        funRef_(iRef_).name = repFunNames_{f_}; %#ok
        funRef_(iRef_).description = desc_; %#ok
    end
end


% -------------------------------------------------------------------------
% nested functions
% -------------------------------------------------------------------------

    function [outFilename, shortDesc_] = parseMFile(inFilename, outDir)
        %parseMFile  parses mfile and writes marked up help content
        % 
        % the m-file inFilename is parsed, the help text extracted and
        % marked up, and the result is written to outFilename in outDir
        % 
        
        [~, fileName, ext] = fileparts(inFilename);
        
        fidIn = fopen(inFilename, 'r');
        oncleanup_in = onCleanup( @() fclose(fidIn) );
        outStr_ = '';
        
        % initialize short description
        shortDesc_ = '';
            
        outFilename = fullfile(outDir, [fileName opts.OutSuffix ext]);
        
        
        % initialize variables to ensure they are in the scope of nested
        % functions
        myLine_ = ''; lenLine_ = 0; currentBlock_ = '';
        
        table = struct('nCol', {}, 'curCol', {}, 'nBlank', {}, ...
            'isOpen', {}, 'style', {});
        
        idxOpen_ = 0;
        
        envIsOpen_ = false;
        
        codeIsOpen_ = false;
        
        div = struct('isOpen', false, 'class', '');
        
        % count braces to enable new paragraphs in tables using semi-colons
        bracesCounter_ = struct('curly', 0, 'normal', 0, 'square', 0);
        
        % first line
        nextLine;
        
        % get function name
        name = regexp(myLine_, ...
            '\s*function(?:[\s\w\[\],]+=\s*|\s*)(?<fn>[\w_]+)', 'names');
        if(~isempty(name))
            funName = name.fn;
        else
            funName = '';
        end
        
        % check whether the file contains a function definition
        if ~isempty(funName)
            try
                assert( strcmp(funName, fileName), ['The functionname ' ...
                    '''%s'' and the filename ''%s'' do not match.'], ...
                    funName, fileName);
            catch me
                if opts.FormatErrors
                    rethrow(me);
                else
                    fprintf(2, ['Format error:\n''%s''\nThe file will ' ...
                        'be skipped!'], me.message);
                    outFilename = '';
                    return
                end
            end
            nextLine;
        % if not, use the filename
        else
            try
                assert( lenLine_ > 0 && strcmp(myLine_(1), '%'), [ ...
                    'The file:\n''%s''\ndoes not seem to be a valid file\n' ...
                    'containing help comments.'], fileName);
            catch me
                if opts.FormatErrors
                    rethrow(me);
                else
                    fprintf(2, ['Format error:\n''%s''\nThe file will ' ...
                        'be skipped!\n'], me.message);
                    outFilename = '';
                    return
                end
            end
            funName = fileName;
        end
        
        if ~strncmp(myLine_, '%', 1)
            fprintf(2, ['The current file:\n''%s''\ndoes not contain ' ...
                'help comments. Will be skipped!\n'], inFilename);
            outFilename = '';
            return;
        end
        
        % look for function name title in the help string
        [title_match, tn] = regexpi(myLine_, ...
            ['^%\s*' funName '\s*-{0,1}\s*(?<desc>[^\s].*)'], ...
            'match', 'names');
        
        % verify the title line
        try
            assert( ~isempty(title_match), ['The functionname in the ' ...
                'title line of the help comment:\n\t''%s''\ndoes not ' ...
                'match the function name ''%s''.'], myLine_, funName);
        catch me
            if opts.FormatErrors
                rethrow(me);
            else
                fprintf(2, ['Format error:\n''%s''\nThe file will ' ...
                    'be skipped!\n'], me.message);
                outFilename = '';
                return
            end
        end
        
        % write the header and following description
        appendString( ['%% ' funName] );
        if ~isempty(tn.desc)
            appendString( ['% ' tn.desc] );
            shortDesc_ = tn.desc;
        else
            warning(['The short description is missing ' ...
                'for ''%s''.'], funName);
        end
        
        % get second line of short description
        nextLine;
            
        % loop the help comment
        while ischar(myLine_) && strncmp(myLine_, '%', 1)
            
            
            % check for skip pattern, i.e. if encountered line is skipped
            skip_ = ~all( arrayfun( @(x)isempty(x{1}), ...
                    regexp(myLine_, opts.SkipPattern) ) );
            if skip_ 
                nextLine;
                continue;
            end
            % check for stop pattern, i.e. if encountered function stops
            % execution
            stop_ = ~all( arrayfun( @(x)isempty(x{1}), ...
                    regexp(myLine_, opts.StopPattern) ) );
            if stop_ 
                break;
            end
            % check for empty comment lines
            if strcmp(myLine_, '%')
                
                if envIsOpen_
                    if div.isOpen
                        finishDiv;
                    end
                    if idxOpen_ > 0 && table(end).isOpen.paragraph
                        finishParagraph;
                    end
                elseif codeIsOpen_
                    appendString( '%%' );
                    codeIsOpen_ = false;
                else
                    % add new line
                    appendString( myLine_ )
                end
                
                
                
                % get next line and continue
                nextLine;
                continue;
            end
            
            % replace tabs by spaces
            tabRep;
            % replace function name, ...
            strRep;
            
            % characterize the content of the current line
            [mtch, cont] = regexp(myLine_, [...
                ... heading
                regExpPat.head ... '^% (?<head>[^\s].*):\s*' ...
                '|' ... or
                ... table ...
                regExpPat.table ... '^%(?<nblank>\s{2,})' ...
                ... '-(?<prop>[^:]+):' ... property
                ... '\s*(?<type>[^;]*);' ... type
                ... '\s*(?<text>[^\s]?.*)$' ... cell content
                '|' ... or
                ... note box
                regExpPat.info ...'^%(?<nblank>\s*)' ...
                ... '//\s*(?<info_title>[^\s][^:]*):?' ... title line
                ... '\s*(?<info_text>[^\s:][^:]*)$' ... content of title line
                '|' ... or
                regExpPat.source ... '^%(?<nblank>\s*)' ...
                ... '>(?<source>.*)$' ... source code
                '|' ... or
                regExpPat.text ...'^%(?<nblank>\s{2,})' ...
                ... '(?<text>[^\s].*)$' ... plain text
                ], 'match', 'names');
            cont.nblank = length(cont.nblank);
            
            if isempty(cont.source)
                
                if codeIsOpen_
                    appendString( '%%' );
                    codeIsOpen_ = false;
                end
            end
            
            if isempty(mtch)
                if opts.WriteUnmatched
                    fprintf(2, ['The line:\n''%s''\n could not be matched' ...
                        '. Will be written unmodified.\n'], myLine_);
                    appendString( myLine_ );
                else
                    me = MException( ...
                        'The line:\n''%s''\n could not be matched.',  ...
                        myLine_);
                    if opts.FormatErrors
                        throw(me);
                    else
                        fprintf(2, ['Format error:\n''%s''\nThe file will ' ...
                            'be skipped!\n'], me.message);
                        outFilename = '';
                        return
                    end
                end
                nextLine;
                continue;
            end
            
            % specify the found content
            if ~isempty(cont.head)
                % current line contains a header
                
                % close all open tables, info boxes, ...
                finishOpen;
                
                appendString( ['%% ' cont.head] );
                
                currentBlock_ = cont.head;
                
                nextLine;
                continue;
            elseif ~isempty(cont.prop)
                % current line is table entry
                
                cont.nblank = (cont.nblank + opts.TableBlankCorrection);
                
                while idxOpen_ > 0 && table(end).nBlank > cont.nblank
                    finishCurrent;
                end
                
                if idxOpen_ > 0 && table(end).nBlank == cont.nblank
                    if table(end).isOpen.row
                        finishRow;
                    end
                else
                    if idxOpen_ > 0
                        startParagraph;
                    end
                    startTable;
                end
                
                prop = formatStr(cont.prop);
                if length(prop) > 1 
                    error('property entries should not contain semi-colons')
                end
                prop = strrep(prop{1}, ',', '<br>');
                
                % start row and write first cell in this row
                startRow;
                startCell;
                appendString( sprintf('%%%*s%s', ...
                    3*opts.NSpace, ' ', prop) );
                finishCell;
                
                % next cell
                startCell;
                
            elseif ~isempty(cont.info_text) || ...
                    ~isempty(cont.info_title)
                % current line contains an info box
                
                while idxOpen_ > 0 && table(end).nBlank > cont.nblank
                    finishCurrent;
                end
                
                % check whether an info box has already been opened
                if ~envIsOpen_
                    startEnv;
                end
                
                if ~div.isOpen || ~strcmp(div.class, 'info')
                    startDiv('info');
                    if ~isnumeric(cont.info_title)
                        cont.info_text = sprintf('<b>%s</b>   %s', ...
                            cont.info_title , cont.info_text);
                    end
                else
                    if ~isnumeric(cont.info_title)
                        cont.info_text = sprintf('%s %s', ...
                            cont.info_title , cont.info_text);
                    end
                end
                        
                % check if info title is present
                info_text = formatStr( cont.info_text );
                for iT_ = 1:length(info_text)
                    info_text__ = sprintf('%%%*s%s', 2*opts.NSpace, ' ', ...
                        info_text{iT_});
                    appendString( info_text__ );
                end
                
                
            elseif ~isnumeric(cont.source)
                % current line contains source code
                % write without percent sign
                
                codeIsOpen_ = true;
                
                txt = regexp(myLineBak_, ...
                    '^%(?<nblank>\s*)>(?<source>.*)$', ... 
                    'names');
                
                appendString( ['  ' txt.source] );
                
            end
            
            if ~isempty(cont.text)
                
                
                switch lower(currentBlock_)
                    case 'syntax'
                        if ~envIsOpen_
                            startEnv;
                        end
                        if div.isOpen && ~strcmp(div.class, 'syntax')
                            finishDiv;
                        end
                        if ~div.isOpen 
                            startDiv('syntax');
                        end
            
                        txt_ = cont.text;
                        txt_ = {strrep(txt_, '|', '')};
                        txt = '';
                        for iT_ = 1:length(txt_) 
                            txt = [txt sprintf('%% %s <br>', txt_{iT_})]; %#ok<AGROW>
                        end
                    
                    case 'see also'
                        
                        spl = regexp(myLineBak_, ...
                            '^%\s+(?<text>[^\s].*)$', 'names');
                        spl_txt = regexp(spl.text, ',', 'split');
                        
                        % replace function names by links
                        curSeeAlsoName_ = strtrim(spl_txt{1});
                        if any(strcmp(curSeeAlsoName_, repFunNames_))
                            txt = ['% <' curSeeAlsoName_ ...
                                opts.OutSuffix '.html |' ...
                                curSeeAlsoName_ '|>'];
                        else
                            txt = ['% <matlab:doc(''' curSeeAlsoName_ ...
                                ''') ' curSeeAlsoName_ '>'];
                        end
                        for iS_ = 2:length(spl_txt)
                            curSeeAlsoName_ = strtrim(spl_txt{iS_});
                            if isempty(curSeeAlsoName_)
                                continue;
                            end
                            if any(strcmp(curSeeAlsoName_, repFunNames_))
                                txt = [txt ', <' curSeeAlsoName_ ...
                                    opts.OutSuffix '.html |' ...
                                    curSeeAlsoName_ '|>']; %#ok
                            else
                                txt = [txt ', <matlab:doc(''' curSeeAlsoName_ ...
                                    ''') ' curSeeAlsoName_ '>']; %#ok<AGROW>
                            end
                        end
                        clear spl_txt iS_ curSeeAlsoName_
                        
                    otherwise
                        
                        while idxOpen_ > 0 && table(end).nBlank > cont.nblank
                            finishCurrent;
                        end
                        
                        if idxOpen_ > 0
                            % replace semi-colons with new paragraphs if not
                            % inside of braces (),[],{} or apostrophes ''
                            % cont.text = insertNewPar(cont.text);
                            cont.text = formatStr(cont.text);
                            blanks = sprintf('%*s', 3*opts.NSpace, ' ');
                            if div.isOpen
                                finishDiv;
                            end
                            if ~table(end).isOpen.paragraph
                                startParagraph;
                            end
                        else
                            if envIsOpen_ && ~div.isOpen
                                finishEnv;
                            end
                            blanks = ' ';
                        end
                        if iscell(cont.text)
                            txt = cell(1,length(cont.text));
                            for iT_ = 1:length(cont.text)
                                if isempty(cont.text{iT_})
                                    txt{iT_} = '';
                                else
                                    txt{iT_} = sprintf(['%%' blanks '%s'], ...
                                        cont.text{iT_});
                                end
                            end
                        else
                            txt = sprintf(['%%' blanks '%s'], cont.text);
                        end
                end
                
                if iscell(txt)
                    for iT_ = 1:length(txt)
                        
                        if isempty(txt{iT_})
                            continue;
                        end
                        
                        if idxOpen_ > 0 && ~table(end).isOpen.paragraph
                            startParagraph;
                        end
                        
                        appendString( txt{iT_} );
                        
                        if iT_ < length(txt) && ...
                                (bracesCounter_.normal + bracesCounter_.square + ...
                                bracesCounter_.curly) == 0
                            if idxOpen_ > 0 && table(end).isOpen.paragraph
                                finishParagraph;
                            end
                        end
                        
                    end
                else
                    appendString( txt );
                end
                
            end
            
            nextLine;
            
        end

        % fclose(fidIn);
        
        if ~isempty(opts.Footer)
            outStr_ = [outStr_ sprintf('\n%s', opts.Footer)];
        end
        
        fidOut = fopen(outFilename, 'w');
        fprintf(fidOut, '%s', outStr_);
        fclose(fidOut);
        
        % --------------------------------------------------------------- %
        % nested functions
        % --------------------------------------------------------------- %
        
        function nextLine
            if fidIn <0
                fprintf('strange')
            end
            myLine_ = fgetl( fidIn );
            if ischar(myLine_)
                myLine_ = strtrim( myLine_ );
            end
            lenLine_ = length( myLine_ );
        end
       
        function tabRep
            % replace tabs by spaces (as defined by the 'NSpace' option)
            myLine_ = strrep( myLine_, sprintf('\t'), ...
                sprintf('%*s', opts.NSpace , ' ') );
        end
            
        function strRep
           
            myLineBak_ = myLine_;
            
            if opts.RepFunNames
                % add markup to all processed function names
                for iRep = 1:length(repFunNames_)
                    myLine_ = strrep( myLine_, repFunNames_{iRep}, ...
                        ['|' repFunNames_{iRep} '|']);
                end
            else
                % add markup only to current function name
                myLine_ = strrep( myLine_, funName, ['|' funName '|']);
            end
            
            if ~strcmp(funName, upper(funName))
                % replace upper case function names by original name
                myLine_ = strrep( myLine_, upper(funName), ['|' funName '|']);
            end
                        
            % add monospace markup for strings
            tokIdx = regexp(myLine_, ...
                '[\s\W\(\{\[,](''[^'']*'')[\s\W\)\}\],]', 'tokenExtents');
            for iTok = 1:length(tokIdx)
                idx_ = tokIdx{iTok} + 2*(iTok - 1);
                myLine_ = [myLine_(1 : idx_(1)-1) '|' ...
                    myLine_(idx_(1) : idx_(2)) '|' ...
                    myLine_(idx_(2)+1 : end)];
            end
        end
        
        function appendString( s )
            outStr_ = [outStr_ sprintf('%s\n', s)];
        end
        
        function finishCurrent
            
            if div.isOpen
                finishDiv;
            end
            
            if table(end).isOpen.paragraph
                finishParagraph;
            end
            
            % finish current active html table
            if idxOpen_ > 0
                finishTable;
            end
            
            if idxOpen_ == 0
                finishEnv;
            end
        end
        
        function finishOpen
            if div.isOpen
                finishDiv;
            end
            if idxOpen_ > 0 && table(end).isOpen.paragraph
                finishParagraph;
            end
            % finish all open html tables
            while idxOpen_ > 0;
                finishCurrent;
            end
            if envIsOpen_
                finishEnv;
            end
        end
        
        
        function out = formatStr( in )
            
            % divide into paragraphs when semi-colon appears not in between
            % braces. braces of previous lines are counted as well. a new
            % paragraph reinitializes the counter
            newPar_ = regexp(in, ';', 'split');
            out = {''};
            curOutIdx_ = 1;
            for iN_ = 1:length(newPar_)
                bracesCounter_.normal = bracesCounter_.normal + ...
                    numel(strfind(newPar_{iN_}, '('));
                bracesCounter_.normal = max( 0, bracesCounter_.normal - ...
                    numel(strfind(newPar_{iN_}, ')')) );
                
                bracesCounter_.square = bracesCounter_.square + ...
                    numel(strfind(newPar_{iN_}, '['));
                bracesCounter_.square = max( 0, bracesCounter_.square - ...
                    numel(strfind(newPar_{iN_}, ']')) );
                
                bracesCounter_.curly = bracesCounter_.curly + ...
                    numel(strfind(newPar_{iN_}, '{'));
                bracesCounter_.curly = max( 0, bracesCounter_.curly - ...
                    numel(strfind(newPar_{iN_}, '}')) );
                
                if (bracesCounter_.normal + bracesCounter_.square + ...
                        bracesCounter_.curly) > 0
                    out{curOutIdx_} = [out{curOutIdx_}, newPar_{iN_}, ';'];
                else
                    out{curOutIdx_} = [out{curOutIdx_}, newPar_{iN_}];
                    curOutIdx_ = curOutIdx_ + 1;
                    out{curOutIdx_} = '';
                end
            end
            if isempty(out{end}) && ( ~strcmp(in(end), ';' ) || ...
                    strcmpi(currentBlock_, 'syntax') )
                out(end) = [];
            end
            
            % replace curly braces
            pat = {'{', '}',};
            switch opts.Format
                case 'html'
                    rep = {'&#123;', '&#125;'};
                case {'pdf', 'latex'}
                    rep = {'\{', '\}'};
            end
            
            for iO_ = 1:length(out)
                for i_ = 1:length(pat)
                    out{iO_} = strrep(out{iO_},  pat{i_}, rep{i_});
                end
            end
            
            
            % replace matlab publish markup with its html/latex equivalent
            pat = {'\|', '\*', '_'};
            switch opts.Format
                case 'html'
                    rep = {'tt', 'b', 'i'};
                    left = '<%s>';
                    right = '</%s>';
                case {'pdf', 'latex'}
                    rep = {'tttext', 'bftext', 'ittext'};
                    left = '\\%s\{';
                    right = '\}';
            end
            
            
            for iO_ = 1:length(out)
                for i_ = 1:length(pat)
                    p_ = pat{i_};
                    out{iO_} = regexprep(out{iO_}, ...
                        ['(?<!\w)' p_ '([^\(\)\[\]\{\}' p_ ']+)' p_ '(?!\w)'], ...
                        sprintf([left '$1' right], rep{i_}, rep{i_})); %#ok<AGROW>
                end
            end
            
        end
        
        function startDiv(classID_)
            div_ = sprintf('%%%*s%s', 2*opts.NSpace, ' ', ...
                ['<div class="' classID_ '">']);
            appendString( div_ );
            
            div.isOpen = true;
            div.class = classID_;
        end
        function finishDiv()
            div_ = sprintf('%%%*s%s', 2*opts.NSpace, ' ', '</div>');
            appendString( div_ );
            
            div.isOpen = false;
            div.class = '';
        end
        
        function startParagraph(classID_)
            if nargin == 0 
                classID_ = 'table';
            end
            if table(end).isOpen.paragraph
                finishParagraph;
            end
            par = sprintf('%%%*s%s', 2*opts.NSpace, ' ', ...
                ['<p class="' classID_ '">']);
            appendString( par );
            
            table(end).isOpen.paragraph = true;
        end
        function finishParagraph
            if div.isOpen
                finishDiv;
            end
            par = sprintf('%%%*s%s', 2*opts.NSpace, ' ', '</p>');
            appendString( par );
            
            table(end).isOpen.paragraph = false;
            
            bracesCounter_.normal = 0;
            bracesCounter_.square = 0;
            bracesCounter_.curly = 0;
        end
        
        function startTable(classID_)
            if ~envIsOpen_
                startEnv;
            end
            
            if nargin == 0;
                if idxOpen_ > 0
                    classID_ = 'inner';
                else
                    classID_ = '';
                end
            end
            
            newTable = getTableStyle(classID_);
            newTable.nBlank = cont.nblank;

            table = [table newTable];
            idxOpen_ = length(table);
            
            switch opts.Format
                case 'html'
                    tbl = sprintf('%% %s\n%% %s', ...
                        ['<table cellspacing="0" cellpadding="4" ' ...
                        'width="' table(end).style.width '" ' ...
                        'border="' table(end).style.border '" ' ...
                        'frame="' table(end).style.frame '" ' ...
                        'rules="' table(end).style.rules '" ' ...
                        'class="' table(end).style.class '">'], ...
                        ['<tbody valign="' table(end).style.vAlign '">']);
                case {'pdf', 'latex'}
                    tbl = '% \begin{tabular}';
            end
            appendString( tbl );
        end
        
        function finishTable
            % close open row
            if table(end).isOpen.row
                finishRow;
            end
            
            switch opts.Format
                case 'html'
                    env = sprintf('%% %s\n%% %s', '</tbody>', '</table>');
                case {'pdf', 'latex'}
                    env = '\end{tabular}';
            end
            appendString( env );
            
            table(end) = [];
            idxOpen_ = length(table);
        end
        
        function startRow
            
            switch opts.Format
                case 'html'
                    el = sprintf('%%%*s%s', opts.NSpace, ' ',  ...
                        ['<tr bgcolor="' table(end).style.rowCol '">']);
                case {'pdf', 'latex'}
                    el = '';
            end
            appendString( el );
            table(end).isOpen.row = true;
        end
        function finishRow
            % write remaining cells
            if table(end).isOpen.cell
                finishCell;
            end
            while table(end).curCol < table(end).nCol
                startCell;
                finishCell;
            end
            
            % place end of row mark
            switch opts.Format
                case 'html'
                    el = sprintf('%%%*s%s', opts.NSpace, ' ', ...
                        '</tr>');
                case {'pdf', 'latex'}
                    el = '\\ \hline';
            end
            appendString( el );
            table(end).isOpen.row = false;
            table(end).curCol = 0;
        end
        
        function startCell
            switch opts.Format
                case 'html'
                    el = sprintf('%%%*s%s', opts.NSpace, ' ', ...
                        '<td>');
                case {'pdf', 'latex'}
                    el = '';
            end
            appendString( el );
            table(end).isOpen.cell = true;
            table(end).curCol = table(end).curCol + 1;
            startParagraph('table');
        end
        function finishCell
            if ~table(end).isOpen.cell
                return;
            end
            if table(end).isOpen.paragraph
                finishParagraph;
            end
            switch opts.Format
                case 'html'
                    el = '</td>';
                case {'pdf', 'latex'}
                    el = '& ';
            end
            appendString( sprintf('%%%*s%s', opts.NSpace, ' ', ...
                         el) );
            table(end).isOpen.cell = false;
        end
        
        function startEnv
            switch opts.Format
                case 'html'
                    env = '<html>';
                case {'pdf', 'latex'}
                    env = '\begin{tabular}';
            end
            appendString( sprintf('%%\n%% %s', env) );
            envIsOpen_ = true;
        end
        
        function finishEnv
            
            switch opts.Format
                case 'html'
                    env = '</html>';
                case {'pdf', 'latex'}
                    env = '\end{tabular}';
            end
            appendString( sprintf('%% %s\n%%', env) );
            envIsOpen_ = false;
        end
        
        function newT = getTableStyle(classID_)

            if nargin == 0
                classID_ = '';
            end

            newT.curCol = 0;
            newT.isOpen = struct('row', false, 'cell', false, 'paragraph', false);
            newT.nBlank = [];
            
            newT.nCol = 2;
            newT.style.rowCol = '#F2F2F2';
            newT.style.borderCol = '#B2B2B2';
            newT.style.border = '1';
            newT.style.width = '';%'100%';
            newT.style.vAlign = 'top';
            newT.style.rules = 'none';
            newT.style.frame = 'box';
            newT.style.class = classID_;
        end
    end


end
