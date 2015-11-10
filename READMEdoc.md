# sssMORdoc
Documentation for sssMOR, a sparse state-space, model order reduction toolbox developed at the Chair of Automatic
Control, Technische Universitaet Muenchen.

*Programmed and tested with: MATLAB R2015b*

_________


### Table of contents

[TOC]


## sssMORdoc on GitLab
On your computer, clone the sssMORdoc repository with SourceTree. The sssMORdoc repository includes the sssMOR
repository as a submodule  and the sssMOR repository includes the sss repository as submodule. That means that when you
clone the sssMORdoc repository you automatically  also clone the sssMOR and the sss repository!

The submodule repositories (sssMOR and sss) will be inside the working directory of the sssMORdoc and you can treat them
as normal repositories: You can change the files in the submodule repositories, and when you are happy with the results,
you can commit and push those changes. You can also pull the newest version of the repositories.

Keep in mind that if you update a submodule (e.g. you pull changes or you make changes to it), the "parent" repository
will see this as a change  and has to be updated by commiting them. Please don't forget push commited changes. :)

The reason why the doc functions are not included in sssMOR is that they should be kept private (many functions are
based on code from Simtech  in Stuttgart). Using this structure, it will be possible to work on the documentation  and
on the toolbox **in the same directory**.

In sssMOR, the **headerTemplate** should be written/commented such that developers of sssMOR know how to format the
header in order to produce the desired outcome.


## How to create the documentation
After having cloned the sssMORdoc repository, add it (and all of its subfolders) to the Matlab Path. You can do this
from the Current Folder Window in Matlab by doing a right-click on the repository folder and selecting 
"Add to the path > Selected Folders and Subfolders". Alternatively, you can do the same from  "Home > Set Path".

After having added sssMORdoc to your path, run "publishHelp" from the Command Window. This will generate automatically
the HTML documentation. Go to the Matlab Documentation (type "doc" in the Command Window) and you should now be able to
see, on the lower right corner (under "Supplemental Software"), the "sssMOR Toolbox" documentation.

The "publishHelp.m" function will update the following files everytime it is run: 
- doc/source/functions_index.m
- doc/source/functions/*.m
- doc/html/helptoc.xml
- doc/html/*.html


## Documentation Workflow
If you want to make changes to the **any** documentation file (e.g. the header of a function) then follow this steps:

1. Open the Matlab Documentation ("doc" in Command Window) and search for the HTML page of the file that you want to
change.

2. Switch to the Command Window and run "publishHelp".

3. Switch back to the Matlab Documentation Window and then press F5 on your Keyboard (if you are on a Mac, you might
have to press Command + R). This will refresh the page in the Matlab Documentation Window. This way you see to
actualized version of the page.

4. Make changes to the file you want to change (e.g. header of a function)

5. Repeat steps 2, 3 and 4 iteratively until you are happy with the results.

6. Commit and push changes to the corresponding repositories ;)

(7. Update, commit and push parent repositories in case the changes where made to a submodule)

This workflow ensures that you are always aware of how the changes you make are going to look when we publish them.

So for example if you edit the header of a function in the sss repository and you are already happy with the results
after having checked the published HTML (steps 2-4), then you have to:

1. Commit and push the changes of the sss function header to the sss repository

2. Commit and push the update changes of the sss submodule to the sssMOR repository

3. Commit and push the update changes of the sssMOR submodule to the sssMORdoc repository


## Release Workflow
ToDo

## Project management
Here is a list of important aspects for the development of the toolbox that should be kept available to everybody

### Core questions/to dos
Here is a list of **central** aspect/questions that have to be clarified asap
- [ ] In what format (zip, MATLAB toolbox, app, ...) is the toolbox going to be available for download on the homepage?
- [ ] What type of license should we use? (compare to tax license)

### Team & Assignments

- Lisa:     sss, sssMOR, unittests
- Jorge:    sss, sssMOR, demos
- Rodrigo:  documentation, release
- Siyang:   documentation, release