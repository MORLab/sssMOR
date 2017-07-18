%% sssMOR App - Sparse State-Space and Model Order Reduction Matlab Application
%
% The *sssMOR App* is a Matlab Application for Model Order Reduction,
% which uses the functions from the *sss Toolbox* and from the *sssMOR Toolbox*.  
%
% The main menu of the user interface contains the five items *About*,
% *Loading and Setting up Models*, *Model Order Reduction*, *Postprocessing
% and Visualization* and *System Analysis*. The sequence of these five items
% follows the usual workflow for Model Order Reduction, which starts with
% loading the desired model, continues with the reduction and ends with an
% evaluation of the results.
%
% In the following, the main functionality of the Application
% is demonstrated with an example.
%
%% Loading and Setting up Models
%
% This tab provides basic functionalities to manage the variables in the
% workspace. State-Space Models can either be loaded directly from a data
% file (using *load sss-model*) or can be composed using matrices from the 
% workspace (using *load matrices*). In our case we are using the button 
% *load sss-model* to load the model named *beam* to the Matlab workspace.
%
% <<img/GUIdoc_1.png>>
%
% Now that we have a State-Space Model we can work with, we want to reduce
% the model order using different reduction algorithms.
%
%% Model Order Reduction
%
% In this tab, classic model reduction methods such as *Balanced Truncation*,
% *Modal Reduction* and *Krylov Subspace Methods* can be used for reducing
% sparse state-space objects.
%
% *Balancing & Truncation*
%
% First we want to reduce our model using balancing and truncation. We
% select q = 50 as order for the reduced model and save the result of the
% reduction under the name *sysr_beam_tbr*.
%
% <<img/GUIdoc_2.png>>
%
% *Modal*
%
% Next, we reduce our model using the modal reduction technique. Since we
% want to compare the results from the different algorithms in the end, we
% choose again q = 50 as order for the reduced model. Further, we select
% the *Smallest Magnitude* Option, specify *sysr_beam_modal* as name for
% the reduced model and finally click the *Reduce* button.
%
% <<img/GUIdoc_3.png>>
%
% *Krylov*
%
% The last reduction technique we want to use is the Rational Krylov
% subspace method. For this method expansion points have to be specified.
% For our example, we simply choose the default values that are already set
% in the table, but change the number of moments in the first row from 2 to
% the value 46, because we want the order of the reduced model to be equal
% to 50 again. The result is saved to the workspace using the name
% *sysr_beam_rk*.
%
% <<img/GUIdoc_4.png>>
%
%% Postprocessing and Visualisation
%
% After we created the three reduced models, we now of course want to know
% how good the reduced models approximates the behavior of the original
% model. Therefore, we first compare the bode diagram of the original and
% the reduced models. To do this, we simply move the four models from the
% list on the left side (all models in the workspace) to the list on the 
% right side (models selected for the graph) and then click the
% *Plot* button.
%
% <<img/GUIdoc_5.png>>
%
%% System Analysis
%
% The bode diagram provided us with a first visual impression on how good
% our reduced models are. But now we want to compare the models in terms of
% real numbers. This can be done under the *System Analysis* menu item. We
% select the original model on the left side, one of the reduced models
% on the right side and then click some of the *Calculate* buttons,
% depending on what kind information we want to get to know about the
% selected models.
%
% <<img/GUIdoc_6.png>>
%
%%
% <html>
%   <hr>
%   <p class="copy">&copy; 2015-2017 RT Technische Universit&auml;t M&uuml;nchen
%        <tt class="minicdot">&#149;</tt>
%        <a href="https://www.rt.mw.tum.de/?sssMOR">Website</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:txts/LICENSE.txt">License</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:txts/RELEASE.txt">Release Notes</a>
%   </p>
% <div>
% <table>
%  <tr>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%             <img src="img/logo_sssMOR_long.png" alt="sssMOR_Logo" height="40px">
%      </td>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%      <img src="img/MORLAB_Logo.jpg" alt="MORLAB_Logo" height="40px"></td>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%      <img src="img/Logo_Textzusatz_rechts_engl_Chair.png" alt="RT_Logo" height="40px"></td>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%      <img src="img/TUM-logo.png" alt="TUM_Logo" height="40px"></td>
%  </tr>
% </table>
% </div>
% </html>
