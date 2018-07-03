function [vessel_t,totElements,noise_t,next_t,other_t,maskVessels]=calculateTimeValues(dataShifted,classA,classB,classC,classD)
%function [vessel_t,totElements,noise_t,next_t,other_t,maskVessels]=calculateTimeValues(dataShifted,classA,classB,classC,classD)
%------ Calculate the time values given the data and the classes of vessels, tissue and boundaries around
%------ the vessels
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


%-------------- Regular size checks -----
[rows,cols,levs,timeSamples]                                = size(dataShifted);                %------- Dimensions of data
RCD                                                         = rows*cols*levs ;                  %------- RCD as total size of elements
[rowsM,colsM,levsM,timeSamplesM]                            = size(classA);                     %------- Dimensions of data

if nargin==2
    data                                                    = reshape(dataShifted,[RCD timeSamples]);

    %----- case when only one mask with different values is given
    %----- the output will be a matrix of size  [timeSamples x classes]
    %----- **** 11 May 2006 **** previous rules causes problems when all classes are expected, if no class is produced use zero output
    ClassAStr                                               = classA(:);
    maxClass                                                = max(ClassAStr);        
    minClass                                                = max(0,min(ClassAStr));
    if timeSamplesM==1
        %----- this is the case when a SINGLE thresholding is performed for the first time sample, if it has more dimensions then a mask
        %----- for each time sample has been provided
        for counterClass                                    = minClass:maxClass
            currentClass                                    = (ClassAStr==counterClass);
            totElements(1,maxClass-counterClass+1)          = sum(currentClass);
            if totElements(1,maxClass-counterClass+1)==0
                vessel_t(:,maxClass-counterClass+1)         = zeros(timeSamples,1);
            else
                vessel_t(:,maxClass-counterClass+1)         = (mean(data(currentClass,:)))';
            end
        end
        totElements                                         = repmat(totElements,[timeSamples 1]);
    else
        mask                                                = reshape(classA,[RCD timeSamples]);
        for counterClass=minClass:maxClass
            currentClass                                    = (mask==counterClass);
            totElements(maxClass-counterClass+1,:)          = (sum(currentClass));
            tempData                                        = sum(data.*currentClass);
            tempElements                                    = totElements(maxClass-counterClass+1,:);tempElements(tempElements==0)=inf;
            vessel_t(:,maxClass-counterClass+1)             = (tempData)./tempElements;
        end
        totElements                                         = totElements';
    end
   
% 
%     maxClass=max(max(max(classA)));        minClass=max(1,min(min(min(classA))));
%         %----- discard ZEROS, that is; if a mask with zeros and ones is received it is considered as having ONE class
%         vessel_t=[]; totElements=[];
%         for counterClass=minClass:maxClass
%             mask_t1=(classA==counterClass);
%             totElements=[sum(sum(sum(mask_t1))) totElements] ;
%             if totElements==0
%                 vessel_t(:,counterClass) = zeros( timeSamples,1);
%             else
%                 %------ the values are appended to the LEFT, that is, highest value of the mask will occupy first column
%                 %------ and the lowest will occupy the last column
%                 mask=repmat(mask_t1,[1 1 1 timeSamples]);
%                 vessel_t =[squeeze(sum(sum(sum(dataShifted.*mask))))/totElements(1)  vessel_t ];
%             end
%         end
% 
% 


else
    %----- case when several masks are received


    %-----------------
    %-----------------Extend the masks now---------------------------------------------------------
    maskVessels         = repmat(classA,[1 1 1 timeSamples]);
    maskNoise           = repmat(classB,[1 1 1 timeSamples]);
    maskNextToVessels   = repmat(classC,[1 1 1 timeSamples]);
    maskOther           = repmat(classD,[1 1 1 timeSamples]);
    totVessels          = squeeze(sum(sum(sum( maskVessels))));
    totNoise            = squeeze(sum(sum(sum( maskNoise))));
    totNext             = squeeze(sum(sum(sum( maskNextToVessels))));
    totOther            = squeeze(sum(sum(sum( maskOther))));

    %------ obtain numerical averages for each time and each class
    if totVessels==0    vessel_t=zeros( timeSamples,1); else  vessel_t  = squeeze(sum(sum(sum(dataShifted.*maskVessels))))./totVessels;   end
    if totNoise==0      noise_t=zeros( timeSamples,1);  else  noise_t   = squeeze(sum(sum(sum(dataShifted.*maskNoise))))./totNoise;        end
    if totNext==0       next_t=zeros( timeSamples,1);   else  next_t    = squeeze(sum(sum(sum(dataShifted.*maskNextToVessels))))./totNext;   end
    if totOther==0      other_t=zeros( timeSamples,1);  else  other_t   = squeeze(sum(sum(sum(dataShifted.*maskOther))))./totOther;         end


end
% It may be possible to modify this function to receive a variable number of classes and then generate a variable number
% of  time values, but for the moment leave as it is