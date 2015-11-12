%% Matlab Markup
% Examples of Matlab Markup
% 

%% MATLAB(R) Code
% 
%   for i = 1:10
%       disp x
%   end
% 

%% Hello World and Plot
% 
%   disp('Hello, world')
% 
%   t = 0:0.1:10;
% 
%   x = sin(t);
% 
%   plot(t,x);
% 

disp('Hello, world')

t = 0:0.1:10;

x = sin(t);

plot(t,x);

%% Numbered List
% 
% # item1
% # *item2*
% # _item3_
% # |item4|
%

%% Bulleted List
% 
% * item1
% * *item2*
% * _item3_
% * |item4|
% 

%% LaTeX
% 
% LaTeX can be inserted inline: $x^2+e^{\pi i}$
% 
% Or it can also be inserted as a block:
%
% $$e^{\pi i} + 1 = 0$$
%
% so that it isn't inline with the text.

%% Trademarks
% 
% TEXT(TM) 
%
% TEXT(R)
% 
% Matlab(R)
% 
% 

%% Preformated Text
% 
% Not Preformated
% TEXT
% 
%  Preformated
%  TEXT
% 

%% Hyperlink to Local Files
%  
% <rt.gif RT Logo>
% 

%% Matlab Hyperlink
% 
% <https://www.rt.mw.tum.de/startseite-rt/ MORLab>

%% Matlab Dynamic Hyperlinks (code)
% 
% <matlab:disp('hola%20que%20tal') Display "hola que tal" in Command Window>
% 
% <matlab:clear%20all,clc,close%20all Clean up Command Window and Workspace>
% 
% <matlab:clear%20all,clc,load('build'),sys=sss(A,B,C),bode(sys)   Test the sss-Bode Function with the building model>
%
% <matlab:load('build') 1st step: Load building model data>
% 
% <matlab:sys=sss(A,B,C) 2nd step: Load model into sss object>
% 
% <matlab:bode(sys) 3rd step: Plot bode diagram of system>

%% Images in "/html" Folder (Local)
% 
% *JPEG*
% 
% <<iss1.jpg>>
% 
% *GIF*
% 
% <<OCISLY.gif>>
% 

%% Image from Relative Path 
% 
% <<../../html/img/gyro.jpg>>
% 
 
%% Image from URL
% 
% <<https://static-s.aa-cdn.net/img/ios/528998089/f5e700c4756ad6319082a95628f603ce?v=1>>
% 
 