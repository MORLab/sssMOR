function [output1,output2] = headerTemplate(input1,varargin)
% HEADERTEMPLATE - sssMOR toolbox header template
% 
% Syntax:
%		output1				= HEADERTEMPLATE(input1)
%		output1				= HEADERTEMPLATE(input1,input2)
%		[output1,output2]	= HEADERTEMPLATE(input1,input2,Opts)
% 
% Description:
%       This header file is to be copied and pasted at the beginning
%       of every script and function within the sssMORLab toolbox for reference.
% 
%       The description should be kept as short as possible but as long as
%       required to give a good impression of what the function/script is 
%       about and how it should be used.
%       
%       You can use Matlab Markup in normal paragraph text, but not inside 
%       HTML tables (for example the tables from input/output arguments art
%       HTML tables)
%
%       Matlab Markup doesn't 
%       recognize line breaks.
%
%		//Note: the building.mat file for the example can be found in the
%		benchmarks folder or you can download them from the SLICOT website
%
%       It only recognizes a new paragraph when theres is an empty line in
%       between.
%
%       Due to the parsing of the headers by the publishHelp function, the
%       "preformatted text" and "code" feature of Matlab Markup can't be
%       used normally. To display code use the ">" as shown below
%
%       To display Matlab code use the ">" marker, for example
%
%>  x = 2*pi*(0:0.01:1)
%>  y = sin(x)
% 
%       This generates inline |MONOSPACED Text|
% 
%       This generates _ITALIC Text_
% 
%       This generates *BOLD Text*
% 
%       This generates a <https://www.rt.mw.tum.de URL>
%
%       This generates inline $LaTeX: x^2+e^{\pi i}$
% 
%       To insert a bulleted or numbered list, leave an empty line before
%       the list:
% 
%       * ITEM1
%       * ITEM2
% 
%       # ITEM1
%       # ITEM2
%
%       To insert an image, leave an empty line and then put the filename
%       between ""<<_>>". Notice that the image can be local or online. If
%       the image is local, it has to be stored inside of doc/html/img. The
%       path to the local image is the path between the generated HTML file
%       and the image.
%
%        *Local image:*
%
%       <<img/iss1.jpg>>
%       
%       *Internet image:*
% 
%       <<https://www.rt.mw.tum.de/fileadmin/w00bhf/www/_migrated_pics/Logo_blau_Homepage_01.gif>>
%
% Input Arguments:
%		*Required Input Arguments:*
%		-input1: 		Description of in1
%		-input2: 		Description of in2
%
%		*Optional Input Arguments:*
%		-input3:   		Description of in3
%						//Note: note inside a table must be on new line
%		-Opts:	 		Description e.g.: A structure containing following fields
%			-.field1:  	description of Opts.field1 function;
% 						[opt1 / opt2 / {default opt}]
%			-.field2:  	description of Opts.field2 function;
% 						[opt1 / {opt2} / opt3]
%			-.field3:  	description of Opts.field3 function;
% 						[{opt1} / opt2 / opt3]
%						or line break on a new line
%						//Note: note inside a table must be on new line
%
% Output Arguments:
%       -output1: Description 
%       -output2: Description
%
% Examples:
%		This code loads a benchmark model included in the toolbox
%		and plots its bode diagram using the sparse state space class:
%
%> load building; 
%> sys = sss(A,B,C);
%> bode(sys);
%
%//Note: the building.mat file for the example can be found in the benchmarks
%		folder or you can download them from the SLICOT website
%
% See Also: 
%		arnoldi
%
% References:
%		* *[1] Grimme (1997)*, Krylov projection methods for model reduction
%		* *[2] Antoulas (2005)*, Approximation of large-scale dynamical systems
%		* *[3] Giraud (2005)*, The loss of orthogonality in the Gram-Schmidt...
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Author 1, Author 2
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------