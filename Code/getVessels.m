function [classMask] = getVessels(dataIn,g_marker,g_noise,border)
%function  [classMask] = getVessels(dataIn,g_marker,g_noise,border)
%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro----------------------------------------------------
%------             Postdoc  Sheffield University    ----------------------------------------------------
%------  13 January 2006   ------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------

%------ Segment into regions according to the intensity values
%------ dataIn is the data, normally a 3D set
%------ input:      g_marker    defines vessel intesity,       it is an upper threshold
%------             g_noise     defines background intensity,  it is a lower threshold
%------             border      is a surrounding region for the vessels it can be one thick border or
%------             two concentric surrounding sections
%------ output      one classMask with several values:
%------             4 - vessel
%------             3 - noise
%------             2 - next to vessel
%------             1 - next to noise
%------             0 - everything else


%------ regular dimension checks
[rows,cols,levs,timeSamples]=size(dataIn);                                          %------- Dimensions of data
RCD=rows*cols*levs ;                                                                %------- RCD as total size of elements


if nargin==3
    border=g_noise;
    varThresholds=g_marker;
end

%-----  define a mask that will be used to eliminate isolated pixels with function PIXVSN.m
if levs>1 pixMask=[3 3 3]; else pixMask=[3 3]; end
%-----  define a Gaussian mask that will be used to smooth the final results
if levs>1
    convKernel=gaussF(border,border,ceil(border/20));
else
    convKernel=gaussF(border,border,1);
end
%------ if Time Samples>1 then a variable threshold in time is possible, it should be one [timex2] variable
%------ the first column is the high threshold, the second is the low threshold

if timeSamples==1
    %-------------------------------------------------------------------------------------------------
    %----- Just **** ONE **** time sample of the data, get the classes accordingly ------------------- 
    %-----  Check if that the marker is higher than the noise
    if g_marker<g_noise
        temp=g_marker;    g_marker=g_noise;    g_noise=temp;
    end

    vesselReg=      (dataIn>g_marker);                                              %----- above a threshold --> vessels with marker
    noiseReg=       (dataIn<g_noise);                                               %----- below a threshold --> background (noise)
    noiseReg=       pixvsn(noiseReg,pixMask,2)-1;                                   %----- fill in the holes erodes a bit the borders

    %------ dilate  vesselReg with a gaussian mask
    borderReg=(convn(vesselReg,convKernel,'same'));                                 %----- Region close to marker --- Gaussian smoothing
    vesselReg=(vesselReg+borderReg>0.8)>0;                                          %----- include part of the smooth in *vesselReg* fills holes and adds smoothness but in vesselreg some eroded parts of borderReg are kept
    %--- for future cases:
    %---    if a boundary of thickness = 1 is required use only a zerocrossing of
    %---    vesselReg -0.5
    borderReg1=((borderReg>0.2)-vesselReg)>0;                                       %----- The *border* will be all the smoothed region but without vessels
    borderReg2=((borderReg>0.09)-vesselReg-borderReg1)>0;                           %----- The *border* will be all the smoothed region but without vessels
    noiseReg=(noiseReg+1-borderReg1 - borderReg2-vesselReg)>0;                      %----- the *noise* must not include the border
    %noiseReg=(noiseReg-vesselReg)>0;                                               %----- nor the vessels
    %ClassC=1-(vesselReg+noiseReg+borderReg);                                       %----- everything else
    %if nargout==3 noiseReg=noiseReg+ClassC; end                                    %----- in case 3 classes are requested, merge noise with other

    classMask=vesselReg*4+noiseReg*3+borderReg1*2+borderReg2;

else
    %-------------------------------------------------------------------------------------------------
    %------ there are several time samples, proceed iteratively to get classes for each of them
    for counterTime=1:timeSamples
        classMask(:,:,:,counterTime)=getVessels(dataIn(:,:,:,counterTime),varThresholds(counterTime,1),varThresholds(counterTime,2),border);
    end
end
    %-------------------------------------------------------------------------------------------------

