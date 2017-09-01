sssMOR
=======

A sparse state-space, model order reduction toolbox developed at the Chair of
Automatic Control, Technische Universität München.

For more information, type `doc` in the command window or visit http://www.rt.mw.tum.de/?sssMOR. Check out also our demo by typing `sssMOR_gettingStarted` in the command window

***
*Programmed with:* MATLAB R2015b

*Tested on:* MATLAB R2014b, R2015b, R2016b (both Windows 7 and Ubuntu 16.04.1 LTS)

*Some functions require:* Control System Toolbox, Optimization Toolbox.

> Note: The Gitlab repositories hosting the sss and sssMOR projects will become public soon, making it even easier for you to contribute to the sssMOR project.
>Sign up for our newsletter under https://lists.lrz.de/mailman/listinfo/sssmor to stay up to date.

***
Copyright
----------
This toolbox is developed by [MORLab](https://www.rt.mw.tum.de/?morlab), the model reduction lab at the [Chair of Automatic Control](https://www.rt.mw.tum.de/en/home/).

***
Acknowledgements
-----------------
The developing team is thankful to all the research assistants and students at [MORLab](https://www.rt.mw.tum.de/?morlab) that have contributed at creating and developing the sssMOR toolbox since 2010.

The team of [Morembs](http://www.itm.uni-stuttgart.de/research/model_reduction/MOREMBS_en.php), a model reduction software for elastic multibody systems, is sincerely acknowledged for the support in the automated generation of the documentation for the toolbox.

***
Developing guidelines
----------------------

We hope that you enjoy the the toolbox and would like to contribute by extending its capability.
To make sure that the developing does not get out of hand, we prepared a few guidelines that we ask you to follow.


### Folder structure
The folder structure of the toolbox is as follows
- **sssMOR** (main folder)
 	- **app**
 	- **demos**
 	- **doc**
	- **src** (source code)
		- **extras**
		- **LyapunovEq** 
		- **MOR** (reduction algorithms)
		    - **@ssRed** (class definition for reduced objects)
		    - **classic**
		    - **stateOfTheArt**
	- **test**


### Documentation
To automatically generated the documentation for the toolbox from the function headers, type `publishDoc('sssMOR')` in the command window. Make sure to format the function headers according to the ``headerTemplate.m`` provided. To publish the documentation for a single function, use syntax `publishFunction('function name')`.
