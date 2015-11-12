%% Matlab Markup with HTML
% Examples of HTML in Matlab Markup
% 

%% HTML Headings
% 
% <html>
% <pre style="padding:10px; border:1px solid #d3d3d3; background:#f7f7f7;">
% <h1>This is a h1 heading</h1>
% <h2>This is a h2 heading </h2>
% <h3>This is a h3 heading</h3>
% <h4>This is a h4 heading</h4>
% <h5>This is a h5 heading</h5>
% <h6>This is a h6 heading</h6>
% </pre>
% </html>
% 

%% HTML Paragraphs
% 
% <html>
% <p>
% <b>This is a HTML paragraph.</b> 
% Lorem ipsum dolor sit amet, consectetur adipiscing elit. Quisque
% vestibulum malesuada ex et ultrices. Etiam metus metus, facilisis quis
% rutrum eget, vulputate a urna. Praesent ipsum augue, ultrices vel mi
% sagittis, pellentesque blandit tellus. Phasellus vitae magna pulvinar,
% laoreet dui eget, vestibulum nibh. Proin in risus et arcu tempor sodales.
% Nullam efficitur, libero id elementum varius, diam odio tempus risus,
% quis porttitor mi est nec nisi. Praesent est risus, consequat vitae
% vestibulum in, faucibus at mi. Donec mattis arcu a auctor elementum.
% Suspendisse vestibulum faucibus gravida. Proin metus arcu, venenatis sed
% nulla at, placerat vehicula felis. Nunc elit massa, bibendum.
% </p>
% <p>
% <b>This is another HTML paragraph.</b>
% Lorem ipsum dolor sit amet, consectetur adipiscing elit.<br> Quisque
% vestibulum malesuada ex et ultrices. Etiam metus metus, facilisis quis
% rutrum eget, vulputate a urna.<br> Praesent ipsum augue, ultrices vel mi
% sagittis, pellentesque blandit tellus.<br> Phasellus vitae magna pulvinar,
% laoreet dui eget, vestibulum nibh. Proin in risus et arcu tempor sodales.
% Nullam efficitur, libero id elementum varius, diam odio tempus risus,
% quis porttitor mi est nec nisi.<br> Praesent est risus, consequat vitae
% vestibulum in, faucibus at mi.<br> Donec mattis arcu a auctor elementum.
% Suspendisse vestibulum faucibus gravida.<br> Proin metus arcu, venenatis sed
% nulla at, placerat vehicula felis.<br> Nunc elit massa, bibendum.
% </p>
% <br>
% </html>
% 

%% HTML Hyperlinks
% 
% <html>
% <p>
% <b>HTML Hyperlink:</b><br>
% <a href="https://www.rt.mw.tum.de/startseite-rt/">This is a weblink in HTML</a>
% </p>
% <p>
% <b>HTML Matlab Dynamic Hyperlink (code):</b><br>
% <a href="matlab:disp('hola que tal')">Display "hola que tal" in Command Window</a>
% </p>
% <br>
% </html>
% 

%% HTML Local Image as Hyperlink
% 
% <html>
% <p>
% <a href="https://www.rt.mw.tum.de/startseite-rt/">
% <img src="rt.gif" alt="https://www.rt.mw.tum.de" width="186" height="183" style="border:3px solid blue">
% </a>
% </p>
% <br>
% </html>
% 

%% HTML Table with Styles
% 
% <html>
% <style>
% table {
%     width:100%;
% }
% table, th, td {
%     border: 1px solid black;
%     border-collapse: collapse;
% }
% th, td {
%     padding: 5px;
%     text-align: left;
% }
% table#t01 tr:nth-child(even) {
%     background-color: #eee;
% }
% table#t01 tr:nth-child(odd) {
%    background-color:#fff;
% }
% table#t01 th	{
%     background-color: black;
%     color: white;
% }
% </style>
% <table>
%   <tr>
%     <th>First Name</th>
%     <th>Last Name</th>		
%     <th>Points</th>
%   </tr>
%   <tr>
%     <td>Jill</td>
%     <td>Smith</td>		
%     <td>50</td>
%   </tr>
%   <tr>
%     <td>Eve</td>
%     <td>Jackson</td>		
%     <td>94</td>
%   </tr>
%   <tr>
%     <td>John</td>
%     <td>Doe</td>		
%     <td>80</td>
%   </tr>
% </table>
% <br>
% <table id="t01">
%   <tr>
%     <th>First Name</th>
%     <th>Last Name</th>		
%     <th>Points</th>
%   </tr>
%   <tr>
%     <td>Jill</td>
%     <td>Smith</td>		
%     <td>50</td>
%   </tr>
%   <tr>
%     <td>Eve</td>
%     <td>Jackson</td>		
%     <td>94</td>
%   </tr>
%   <tr>
%     <td>John</td>
%     <td>Doe</td>		
%     <td>80</td>
%   </tr>
% </table>
% <br><br><br><br>
% </html>
% 

%% HTML with Scripts
% 
% <html>
% <style>
% table, td {
%     border: 1px solid black;
% }
% </style>
% <h3>A demonstration of how to access a TABLE element</h3>
% <table id="myTable">
%   <tr>
%     <td>cell 1</td>
%     <td>cell 2</td>
%   </tr>
%   <tr>
%     <td>cell 3</td>
%     <td>cell 4</td>
%   </tr>
% </table>
% <p>Click the button to remove the first row in the table.</p>
% <button onclick="myFunction()">Try it</button>
% <script>
% function myFunction() {
%     var x = document.getElementById("myTable");
%     x.deleteRow(0);
% }
% </script>
% <br><br><br><br>
% </html>
% 

