function release_generator
% RELEASE_GENERATOR - a toolbox release folder generator
% 
% Description:
%       This function copies all relevant toolbox files in a release folder
%       that can then be packaged as a Matlab Toolbox. The release folder
%       is generated in the parent folder of the release_generator
%       function.
% 
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Rodrigo Mancilla
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  07 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% configuration
product_name = 'sssMOR';
release_version = '1.2';
ver_str = strrep(release_version, '.', '_');
out_dir_name = [product_name '_Toolbox' '_release_' ver_str];

%% re-publish doc
fprintf(2,'Do you want to re-publish the documentation?\n[Y]/N: ');
inp = input('','s');
inp = strtrim(inp);
inp = lower(inp);
if isempty(inp)
    inp = 'y';
end
if inp == 'y'
    publishHelp;
    fprintf(2,'Re-publishing documentation\n');
elseif inp == 'n'
    fprintf(2,'Keeping old documentation\n');
else
    fprintf(2,'Wrong input. Aborting function without re-publishing\n');
end

%% generate source and destination directories paths
[proj_dir,~,~] = fileparts(mfilename('fullpath'));
[proj_dir_parent,~,~] = fileparts(proj_dir);
out_dir = [proj_dir_parent filesep out_dir_name];

if exist(out_dir,'dir') == 7
    fprintf(2,'WARNING! Release Package already exists!\nOverwrite? [Y]/N: ');
    inp = input('','s');
    inp = strtrim(inp);
    inp = lower(inp);
    if isempty(inp)
        inp = 'y';
    end
    if inp == 'y'
        rmdir(out_dir,'s');
        fprintf(2,'Overwritting old Release Package\n');
    elseif inp == 'n'
        fprintf(2,'Aborting function without overwritting\n');
        return;
    else
        fprintf(2,'Wrong input. Aborting function without overwritting\n');
    end
end


%% doc
src_root_doc = [proj_dir filesep];

dir_html = ['doc' filesep 'html'];
dir_info_xml = ['doc' filesep 'info.xml'];

doc_cell = {dir_html , dir_info_xml};

%% sssMOR
src_root_sssMOR = [proj_dir filesep 'sssMOR' filesep];

dir_benchmarks = ['benchmarks'];
dir_demos =  ['demos'];
dir_test = ['test'];
dir_MOR_readme = ['README.md'];
dir_MOR_extras = ['software' filesep 'extras'];
dir_MOR = ['software' filesep 'MOR'];
dir_sss_readme = ['software' filesep 'sss' filesep 'README.md'];
dir_sss_class = ['software' filesep 'sss' filesep 'src' filesep '@sss'];
dir_sss_extras = ['software' filesep 'sss' filesep 'src' filesep 'extras'];

sssMOR_cell = {dir_benchmarks , dir_demos , dir_test , dir_MOR_readme...
    dir_MOR_extras , dir_MOR , dir_sss_readme , dir_sss_class , dir_sss_extras};

%% source
src_doc = strcat(src_root_doc,doc_cell);
src_sssMOR = strcat(src_root_sssMOR,sssMOR_cell);

source = [src_doc src_sssMOR];

%% destination
dest_root = [out_dir filesep];
destination = [doc_cell sssMOR_cell];

destination = strcat(dest_root,destination);

%% generate directories
for nCell = 1:length(source)
    nParentDir = fileparts(destination{nCell});
    if exist(nParentDir,'dir') ~= 7
        mkdir(nParentDir);
    end
end

%% copy files
for nCell = 1:length(source)
    copyfile(source{nCell},destination{nCell},'f');
end

%% go to release folder
cd(out_dir);
