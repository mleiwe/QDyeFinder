function [SampleN]=mnl_DetermineSampleSize(Zpercent,N,e)
% Function to calculate the sample size of a given population
% Inputs
% Zpercent - the desired confidence value (e.g. 95% -will correct for
% two-sided calculations
% N - the size of the population
% e - margin of error in percent

%% Pre-calculations
zP=(1-(Zpercent/100))/2;
zVal=norminv(1-zP);
E=e/100;
%% Top Half
Top1=(zVal*zVal)*(0.5*(1-0.5));
Top2=E*E;

Top=Top1/Top2;
%% Bottom Half
Bottom=1+(Top1/(Top2*N));
%% Final Equation
SampleN=Top/Bottom;
