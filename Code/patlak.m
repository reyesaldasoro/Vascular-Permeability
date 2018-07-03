function patlak=patlak(measurementsInTime)

%function patlak=patlak(measurementsInTime)
%------ Calculate the Patlak Parameters
%-------------------------------------------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro                       ----------------------------
%------             University of Sussex                                    ----------------------------
%------             http://www.sussex.ac.uk/Users/ccr22/                    ----------------------------
%------  July 2012                                         ---------------------------------------------
%-------------------------------------------------------------------------------------------------------


%-------------------------------------------------------------------------------------------------------
% This m-file is part of a package to analyse vascular Permeability. For a comprehensive 
% user manual, please visit:
%
%           http://www.sussex.ac.uk/Users/ccr22/permeabilityD.php?idMem=software
%
% Please feel welcome to use, adapt or modify the files. If you can improve
% the performance of any other algorithm please contact me so that we can
% update the packate accordingly.                        
%-------------------------------------------------------------------------------------------------------

%------------ D I S C L A I M E R ----------------------------------------------------------------------
% The m-files, data and information is provided for educational and informational purposes only, 
% and is not intended for professional purposes. The authors shall not be liable for any errors or 
% responsibility for the accuracy, completeness, or usefulness of any information, or method in the 
% content, or for any actions taken in reliance thereon. If you find any errors or problems with the 
% programs we should try to fix them but no promise is made.
%-------------------------------------------------------------------------------------------------------


[t,m]                       = size(measurementsInTime);
%adjustingFunc              = @(k1,k2,k3,t) (k1+k2*t).(1-exp(-k3*t));
ct                          = measurementsInTime(:,2);
cp                          = measurementsInTime(:,1);
intCp                       = cumsum(cp);

for counterTime=1:t
    patlak(counterTime,1)   = (intCp(counterTime))/(cp(counterTime));
end

patlak(:,2)                 = ct./cp;

[q2,iq2]                    = sort(patlak(:,1));
q2(:,2)                     = patlak(iq2,2);
patlak                      = q2;
startpoint                  = [ patlak(1,2) (patlak(t,2)-patlak(1,2))/(patlak(t,1)-patlak(1,1)) 1];
[parameters,model]          = fitCurve(patlak,startpoint);
%parameters
[sse,fittedCurv]            = model(parameters);%sse
%[sse,fittedCurv]=model(startpoint);%sse

patlak(:,3)                 = fittedCurv;
patlak(1:3,4)               = parameters';
