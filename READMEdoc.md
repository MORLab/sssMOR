# sssMORdoc
Documentation for sssMOR, a sparse state-space, model order reduction toolbox developed at the Chair of Automatic Control, Technische Universitaet Muenchen.

*Programmed and tested with: MATLAB R2015b*

___

**Table of contents**

[TOC]

## sssMORdoc on GitLab
On your computer, the files in this repository should be in a **branch** of sssMOR. Include this repository as a **remote** and make sure to push changes to this remote (which you may call "doc").

In this way, you will have **only one** folder on your computer named "sssMOR", and depending what branch you are in you will see the documentation as as folder "doc" or not.

The reason why the doc functions are not included in sssMOR is that they should be kept private (many functions are based on code from Simtech in Stuttgart). Using this branch, it will be possible to work on the documentation **in the same directory** and generate the toolbox to be put for download on the homepage.

In sssMOR, the **headerTemplate** should be written/commented such that developers of sssMOR know how to format the header in order to produce the desired outcome.

## How to create the documentation
In order to automatically generate the HTML documentation, go to the folder "doc" and run the publishHelp.m function.

This will update the following files: 
- doc/source/functions_index.m
- doc/source/functions/*.m
- doc/html/helptoc.xml
- doc/html/*.html

After running the publishHelp.m function, type "doc" in the Matlab Command Window (it is important to stay in the "doc" folder while doing this, because the "info.xml" file is also in that folder)

You should now be able to see on the lower right corner (under "Supplemental Software") the "sssMOR Toolbox" documentation.



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