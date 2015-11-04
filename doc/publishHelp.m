function publishHelp
%publishHelp  generate hmtl-help from function help comments and marked-up
%             help files
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
[002]   06/2013   N. Walz   
---------------------------------------------------------------------------
%}
% END ChangeLog
%% Function Body

productName_ = 'sssMOR';

%% setup directories

% directories containing m-files with help comments

% doc directory
[doc_dir,~,~] = fileparts(mfilename('fullpath'));

% project directory
[proj_dir,~,~] = fileparts(doc_dir);

% (m-files can also be specified directly)
fun_files = { [proj_dir filesep 'software' filesep 'MOR' filesep 'ltiMor' filesep 'classic'] };


% html directory, where all generated html-files are saved to
html_dir = fullfile(doc_dir, 'html', filesep);

if exist( html_dir, 'dir' ) ~= 7
    mkdir( html_dir )
else
    % delete existing html files
    delete( fullfile(html_dir, '*.html') );
    % delete existing png files
    delete( fullfile(html_dir, '*.png') )
end

% source directory, where all marked-up m-files are located
source_dir = fullfile(doc_dir, 'source', filesep);

if exist( source_dir, 'dir' ) ~= 7
    mkdir( source_dir )
    mkdir( fullfile(source_dir, 'functions') )
else
    % delete existing function m-files
    delete( fullfile(source_dir, 'functions', '*.m') );
end


%% Extract help comments from function files

% create custom footer line for html-pages
year_ = datestr(now, 'yyyy');
footer_ = sprintf('%s\n', ...
    '%%', ...
    '% <html>', ...
    '%   <hr>', ...
    [ '%   <p class="copy">&copy; ' year_ ' RT Technische Universität München'], ...
    '%        <tt class="minicdot">&#149;</tt>', ...
    '%        <a href="http://www.rt.mw.tum.de">Website</a>', ...
    '%        <tt class="minicdot">&#149;</tt>', ...
    '%        <a href="file:./LICENSE.txt">Terms of Use</a>', ...
    '%        <tt class="minicdot">&#149;</tt>', ...
    '%        <a href="file:./README.txt">Read me</a>', ...
    '%   </p>', ...
    '% </html>' );

% options for funHelp2MarkUp function
% RODRIGO: options are passed as a struct named 'opts' to funHelp2Markup
opts.OutputDir = fullfile(source_dir, 'functions');
opts.Footer = footer_;

% find all files that shall be marked up
% RODRIGO: create opts.Files field in the 'opts' struct if it doesn't exist
if ~isfield(opts, 'Files') || ~iscell(opts.Files)
    opts.Files = {};
end

% RODRIGO: find all files that shall be marked up by the funHelp2Markup function
for iF_ = 1:length(fun_files)
    fun_ = fun_files{iF_};
    switch exist( fun_ ) %#ok
        case 2  % RODRIGO: if fun_ is an m-file
            s_ = which( fun_ );  % RODRIGO: 'which' displays the full path of m-file 
            opts.Files = [opts.Files; s_];
        case 7  % RODRIGO: if fun_ is a folder
            s_ = what( fun_ ); % RODRIGO: 'what' lists all Matlab relevant files in folder (not subfolders) and returns a struct
            % Ensure a path divider at the end
            if(~strcmp(s_.path(end),filesep))
                s_.path(end+1) = filesep;
            end %if
            for iS_ = 1:numel(s_.m) % RODRIGO: s_.m is the array in the struct with m-files
                opts.Files = [opts.Files; [s_.path s_.m{iS_}]]; % RODRIGO: read all m-files in folder
            end
        otherwise
            fprintf(2, 'Could not find file/directory of:\n\t''%s''\n', ...
                fun_);
            fprintf(2, 'Will be skipped!\n');
    end
end

% extract and markup the help comments of the specified m-files
funRef_ = funHelp2MarkUp( opts );
% RODRIGO: pass desired configurations and files to function funHelp2Markup
% RODRIGO: funHelp2Markup returns filenames of the created 
% and marked-up help m-files including the output directory.


%% Generate helptoc and updated help index files

% go to source directory and get m-files
mDirs_ = {'', 'root', 'sssMOR'; ...
    'getting_started', 'getting_started', 'Getting Started'; ...
    'examples', 'examples', 'Examples'};

helpSource_.functions = struct(...
    'path', fullfile(source_dir, 'functions'), ...
    'item', funRef_, ...
    'entry', 'Functions Reference');

% collect data on files to be published
for d_ = 1:size(mDirs_,1)
    
    s = what( fullfile(source_dir, mDirs_{d_, 1}) );
    
    fld_ = mDirs_{d_, 2};
    ent_ = mDirs_{d_, 3};
    helpSource_.(fld_) = struct('path', s.path, 'item', [], 'entry', ent_);
    
    for f_ = 1:length(s.m)
        
        helpSource_.(fld_).item(f_).file = s.m{f_};
        helpSource_.(fld_).item(f_).name = '';
        
    end
    
end

updateIndexFiles( helpSource_, productName_, footer_, doc_dir );

%% Generate the html-files from the marked-up help m-files



fprintf('Generating the html-helpfiles.\n');
fprintf('\tSource directory:\n\t\t%s\n', doc_dir);

fprintf('\tOutput directory:\n\t\t%s\n', html_dir);
fprintf('Start publishing.\n\n');


% options to passed to the publish function
options.format = 'html';
options.outputDir = html_dir;
options.maxHeight = 800;
options.maxWidth = 600;
options.evalCode = false;
options.stylesheet = [doc_dir filesep 'myStyleSheet.xsl'];


fields_ = fieldnames( helpSource_ );
for f_ = 1:length( fields_ )
    
    curSource_ = helpSource_.(fields_{f_});
    
    % change to the source directory
    % this is necessary in order to be able to execute code if so wanted
    % alternatively, the path could be added to the search path
    % cd( curSource_.path ) % RODRIGO: code in path, no need for cd

    
    for iF_ = 1:length( curSource_.item )
        
        pf_ = publish( curSource_.item(iF_).file, options );
        
        [~, filename, ext] = fileparts( pf_ );
        fprintf('Published ''%s'' to:\n\t%s\n', ...
            curSource_.item(iF_).file, [filename ext]);
    end
end


fprintf('\nFinished publishing.\n\n')

% cd( doc_dir )

%% Build the docsearch database files
% necessary to enable searching in the help browser for the created help
% files
fprintf('Building the search database ...')
builddocsearchdb( html_dir );
fprintf(' done.\n')

