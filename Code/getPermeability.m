function data=getPermeability(data)
%function data=getPermeability(data)
%------ With the masks from the segmentation and the data, calculate
%------             the average value of the class in time
%------             the patlak adjusted values
%------             the regional and time leaks
%------             surface / volume / permeability
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




[data.measurementsInTime,data.elementsPerClass] = calculateTimeValues(data.dataIn,1+data.classA);
if size(data.measurementsInTime,2)==4
    data.measurementsInTime(1,5) =0;
    data.elementsPerClass(1,5) =0;
end

data.patlak             = patlak(data.measurementsInTime(:,1:2));
data.k2                 = data.patlak(2,4);                                 %--- Ki with units in ml/min/ml
data.PSV                = 1e4*data.k2/60;                                   %--- Ki = PS/V with units in 1/s
data.vesVol             = pixvsn(data.classA(:,:,:,1)==4);                  %----- USE the segmented data with one smoothing

%surf_Vox                = 5.2*9.09;                                         %----- this values are determined by the size of the voxel
%vol_Vox                 = 5.2*9.09*5.2;                                     %----- 2.6 for 512 and 5.2x5.2x9.09 for 256x256x6

surf_Vox                = (0.5*data.dimVox_R+0.5*data.dimVox_C)*data.dimVox_Z;                       %----- this values are determined by the size of the voxel
vol_Vox                 = data.dimVox_R*data.dimVox_C*data.dimVox_Z;         %----- 



data.vessVolume         = sum(sum(sum(((data.vesVol)))));
data.vessSurface        = sum(sum(sum(zerocross((data.vesVol)-0.5))));
%data.surf_vol          = ((data.vessSurface*surf_Vox/1e8)/(data.vessVolume*vol_Vox/1e12));     %----- this considers as volume the vessel volume
data.surf_vol           = ((data.vessSurface*surf_Vox/1e8)/(data.RCD*vol_Vox/1e12));            %----- this considers as volume the TOTAL volume
data.permeability       = (1e7*data.k2/60/data.surf_vol);
