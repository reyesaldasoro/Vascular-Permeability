function  dataOut = PermeabilityAnalysis(dataIn)
%function  dataOut = PermeabilityAnalysis(dataIn);
%
%        Main plotting function for the Permeability Data Analysis, this m-file launches the Graphical
%        User Interface through which all the calculations and displays are performed.
%
%-------------------------------------------------------------------------------------------------------
%------- VARARGIN   :   data sets that will be later analysed for vascular permeability the data can be: 
%-------                (1) a series of 3D images with different formats inside one folder 
%-------                (2) a series of folders with 2D images inside another folder 
%-------                (3) a matlab file or series of file
%-------                Besides the usual image formats like TIFF, JPEG, BMP, etc, this file can read
%-------                the Bio-Rad PIC format sometimes associated with multiphoton microscopes
%-------                The data is read with readPermeability data if other than a Matlab file
%------- ARGOUT     :   the data as a 4D matrix with ALL the other parameters
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
      
if nargin == 0        %----- no data received, call readPermeabilityData
    %warning('No data or figure received');
    [dataIn,fileName]                               = readPermeabilityData();
    if isempty(dataIn)
        dataIn=[]; 
        return;
    end
    %createVesselFigure(data);
    %return;
end

%------ Data received may be STRUCT / MATRIX / STRING
%------ if a string is received, it can be for several cases:"
%------         updating the plot       'u_4'
%------         update new parameters   'u_*'
%------         load new data           'datName'
if ~ischar(dataIn)
    if isstruct(dataIn)
        
        if ~isfield(dataIn,'dimVox_R')
            %verify that Dimensions of the Voxel are present
            dataIn.dimVox_R                                           = 1; %[um]
            dataIn.dimVox_C                                           = 1; %[um]
            dataIn.dimVox_Z                                           = 1; %[um]
        end       
        data                                        = createVesselFigure(dataIn);
    else
        % the data is not  a string nor a struct, i.e. a matlab 4D matrix, start calculations
        if ~exist('fileName','var')
            fileName                                = 'data set A';
        end
        data                                        = calculateVesselData(dataIn,fileName);

    end
else
    if dataIn(1:2)~= 'u_'
        %------- start from the beginning, read the data from a name --------------------------------
        
        fileName                                    = dataIn;
        dataIn                                      = readPermeabilityData(dataIn);
        
        if isempty(dataIn)
            
            dataOut=[];
            return;
        end
        if isstruct(dataIn)
            dataIn.datName = fileName;
            PermeabilityAnalysis(dataIn);           
        else    
            calculateVesselData(dataIn,fileName);
        end
    else
        %--------------- - - - - - - - U P D A T E - - - - - - - - - - - - - - - - ------------------
        %-------------------we are in the updating part of the GUI
        %------- cases ----- Prompted by ------------------ Action --------
        %-------   1   -  Cross line is modified        -   move line in data, mask and line plot
        %-------   2   -  Slice or Time are modified    -   re plot slice, mask and line
        %-------   3   -  Erase ROI                     -   capture a ROI and mask it out of processing
        %-------   4   -  PLOT   recalculate            -   whole replot of everything RECALCULATE data
        %-------   5   -  De-select ROI                 -   erase ROI selected in 3
        %-------   6   -  Zoom ROI                      -   discard an outer region and keep only central part
        %-------   7   -  Intensity / patlak plots      -   determine which of all regions predefined should be plotted
        %-------   8   -  mask / time / proj / leak     -   change the mask data to leak, time slice ...
        %-------   9   -  time samples removed          -   remove a couple of samples from the time series ...
        %-------   10  -  ? ? ? ? ? ? ? ?
        %-------   A   -  Remove slices                 -   discard a series of slices, read from new window
        %-------   B   -  Zoom ROI  from coords         -   discard an outer region and keep only central part Read from window

        fig                                         = gcf;
        %------  Get the values that were stored under 'userdata' in the pushbutton for plotting 'pbgraf'
        data                                        = get(fig,'userdata');
        
        
        
        
        %------- Verify that all the fields have been generated
        if ~isfield(data,'PSV')
            data.PSV                                 =  1e4*data.k2/60;
        end
        
        if ~isfield(data,'elementsPerClass')
            [tempVals,data.elementsPerClass]         =  calculateTimeValues(data.dataIn,data.classA);
            clear tempVals;
            set(fig,'userdata',data);
        end
        
        if ~isfield(data,'timeThresholds')
            data.timeThresholds                      =  repmat([data.g_markerA data.g_tissueA],[data.timeSamples 1]);
            set(fig,'userdata',data);
        end
        
        
        

        %------- Get all existing objects
        h = allchild(fig);
        pbgraf       =        findobj(h,'tag','pbgraf');
        texBorder =           findobj(h,'tag','texBorder');
        upThresHandle1 =      findobj(h,'tag','upThresHandle1');
        upThresHandle2 =      findobj(h,'tag','upThresHandle2');
        lowThresHandle1 =     findobj(h,'tag','lowThresHandle1');
        lowThresHandle2 =     findobj(h,'tag','lowThresHandle2');
        plotHandle =          findobj(h,'Position', [0.03    0.83    0.4    0.15]);
        plotMontages =        findobj(h,'Position', [0.53    0.83    0.45    0.15]);
        axisTimeLines =       findobj(h,'Position', [0.53    0.13    0.45    0.16]);
        dataHandle =          findobj(h,'tag','dataHandle');
        classesHandle =       findobj(h,'tag','classesHandle');
        texMarker =           findobj(h,'tag','texMarker');
        texTissue =           findobj(h,'tag','texTissue');
        texMarkerA =          findobj(h,'tag','texMarkerA');
        texTissueA =          findobj(h,'tag','texTissueA');
        lin1  =               findobj(h,'tag','lin1');
        lin2  =               findobj(h,'tag','lin2');
        lin3  =               findobj(h,'tag','lin3');
        lin4  =               findobj(h,'tag','lin4');
        lin5  =               findobj(h,'tag','lin5');
        linRatio  =           findobj(h,'tag','linRatio');
        linElements  =        findobj(h,'tag','linElements');
        linThresholds  =      findobj(h,'tag','linThresholds');
        linPatlak  =          findobj(h,'tag','linPatlak');
        kAText  =             findobj(h,'tag','kAText');
        kBText  =             findobj(h,'tag','kBText');
        k1Text  =             findobj(h,'tag','k1Text');
        k2Text  =             findobj(h,'tag','k2Text');
        k3Text  =             findobj(h,'tag','k3Text');
        permeabText =         findobj(h,'tag','permeabText');
        finVessText  =        findobj(h,'tag','finVessText');
        finTissText  =        findobj(h,'tag','finTissText');
        finRatioText  =       findobj(h,'tag','finRatioText');
        axTim  =              findobj(h,'tag','axisTimeLines');
        texSlice =            findobj(h,'tag','texSlice');
        texTime =             findobj(h,'tag','texTime');
        texXLine  =           findobj(h,'tag','texXLine');
        XLine  =              findobj(h,'tag','XLine');
        XLine2  =             findobj(h,'tag','XLine2');
        XLineP  =             findobj(h,'tag','XLineP');
        XLineP2  =            findobj(h,'tag','XLineP2');
        t_vs_pat =            findobj(h,'tag','t_vs_pat');
        typeOfPlot =          get(t_vs_pat,'value');
        dataToPlot =          get(dataHandle,'cdata');
        sliceToPlot  =        floor(str2num(get(texSlice,'string')));
        colToPlot  =          floor(str2num(get(texXLine,'string')));
        timeToPlot  =         floor(str2num(get(texTime,'string')));
        classesAxis  =        findobj(h,'tag','classesAxis');
        mask_vs_slice =       findobj(h,'tag','mask_vs_slice');
        hist_vs_montage =     findobj(h,'tag','hist_vs_montage');
        typeOfSurf =          get(mask_vs_slice,'value');
        typeOfSurf2 =         get(hist_vs_montage,'value');
        
        

        
        %radButtonAdaptHist  = findobj(gcf,'tag','radButtonAdaptHist'); 
        %adaptiveThres =       get(radButtonAdaptHist,'value');
        %radButtonMulti  =     findobj(gcf,'tag','radButtonMulti'); 
        %multipleTimeThres =   get(radButtonMulti,'value');
        
        % verify that all values are within limits
        timeToPlotN             = min(max(timeToPlot,1),data.timeSamples);
        sliceToPlotN            = min(max(sliceToPlot,1),data.levs);
        colToPlotN              = min(max(colToPlot,1),data.cols);
        
        if timeToPlotN~=timeToPlot
            set(texTime,'string',num2str(timeToPlotN));
            timeToPlot=timeToPlotN;
        end


        if sliceToPlotN~=sliceToPlot
            set(texSlice,'string',num2str(sliceToPlotN));        
            sliceToPlot=sliceToPlotN;
        end

        if colToPlotN~=colToPlot
            set(texXLine,'string',num2str(colToPlotN));        
            colToPlot=colToPlotN;
        end

        
        

        if dataIn(3) == '1'
            %---------------------------------------------------------------------------------------------------------
            %-------------In case the Cross line is modified, update plot and line  ----------------------------------
            %---------------------------------------------------------------------------------------------------------
            if isempty(colToPlot) 
                colToPlot =1;
                set(texXLine,'string',num2str(colToPlot));
            end  
            set(XLine,'xdata',[colToPlot colToPlot]);
            set(XLine2,'xdata',[colToPlot colToPlot]);
            set(XLineP,'ydata',data.dataIn(:,colToPlot,sliceToPlot,timeToPlot));
%             try(data.dataHist);
%                 data.histEqData = data.dataHist(:,:,sliceToPlot,timeToPlot);
%             catch
%                 data.histEqData = (max(max(data.dataIn(:,:,sliceToPlot,timeToPlot))))*adapthisteq(data.dataIn(:,:,sliceToPlot,timeToPlot)/(max(max(data.dataIn(:,:,sliceToPlot,timeToPlot)))));
%             end
%             set(fig,'userdata',data);
%             if (adaptiveThres == 1)|(typeOfSurf == 13)
% 
%                 set(XLineP2,'ydata',data.histEqData(:,colToPlot),'visible','on');
%                 set(plotHandle,'ylim',[min(data.dataIn(:,colToPlot,sliceToPlot,timeToPlot)) 1.05*max(data.histEqData(:,colToPlot))]);
% 
%             else
                set(plotHandle,'ylim',[min(data.dataIn(:,colToPlot,sliceToPlot,timeToPlot)) 1.05*max(data.dataIn(:,colToPlot,sliceToPlot,timeToPlot))]);
                if ~isempty(XLineP2)
                    set(XLineP2,'visible','off');
                end
%             end
            if (typeOfSurf == 2)|(typeOfSurf == 3)|(typeOfSurf == 4)
                clear;gcf;PermeabilityAnalysis('u_8');
            end
            clear;gcf;
            %------  END of UPDATE ** u_1 **  where CROSS LINE is modified -------------------------------------------
        elseif dataIn(3) == '2'
            %---------------------------------------------------------------------------------------------------------
            %----- in case the Slice or time for the data are modified, update surf of data and mask as well as the plot
            %---------------------------------------------------------------------------------------------------------
            %------ the classes
            %set(classesHandle,  'cdata',(data.classA(:,:,sliceToPlot)*4+data.classB(:,:,sliceToPlot)*3+data.classC(:,:,sliceToPlot)*2+data.classD(:,:,sliceToPlot)));
            %            set(classesHandle,  'cdata',data.classA(:,:,sliceToPlot));
            
            if isempty(sliceToPlot)
                sliceToPlot =1;
                set(texSlice,'string',num2str(sliceToPlot));
            end  
            if isempty(timeToPlot)
                timeToPlot =1;
                set(texTime,'string',num2str(timeToPlot));
            end  
                        
            set(dataHandle,     'cdata',data.dataIn(:,:,sliceToPlot,timeToPlot));
            set(XLineP,'ydata',data.dataIn(:,colToPlot,sliceToPlot,timeToPlot));
            set(plotHandle,'ylim',[min(data.dataIn(:,colToPlot,sliceToPlot,timeToPlot)) 1.05*max(data.dataIn(:,colToPlot,sliceToPlot,timeToPlot))]);
            set(fig,'userdata',data);
            if   typeOfSurf>= 1      %(typeOfSurf == 2)|(typeOfSurf == 3)|(typeOfSurf == 4)|(typeOfSurf == 5) |(typeOfSurf == 1)
                clear;gcf;PermeabilityAnalysis('u_8');
            end
            clear;gcf;
            %------  END of UPDATE ** u_2 **  where SLICE of TIME is modified -------------------------------------------
        elseif dataIn(3) == '3'
            %---------------------------------------------------------------------------------------------------------
            %------ this is to capture a region of interest to use as a mask
            %---------------------------------------------------------------------------------------------------------
            classD = (data.classA == 4);
            figure(22);
            maskD = roipoly(dataToPlot);
            close(figure(22))
            maskD = repmat(maskD,[1 1 data.levs]);
            %------- keep the mask as maskD in the data and then use it later on for the correction
            data.maskD = data.maskD|maskD;
            %[vessel_t,tissue_t,next_t,other_t] =       calculateTimeValues(data.dataIn,data.classA,data.classB,data.classC,data.classD);
            %data.measurementsInTime =                      [vessel_t tissue_t next_t other_t];
            set(fig,'userdata',data);
            clear;
            %PermeabilityAnalysis(data);
            %------  END of UPDATE ** u_3 **  where a ROI is selected -------------------------------------------
        elseif dataIn(3) == '4'
            %---------------------------------------------------------------------------------------------------------
            %------- This is the RE PLOTTING  case, recalculate the masks when marker levels are modified
            %------ Update the values with those from the GUI, the ones with blue background [0 1 1] are the
            %------ last ones to be modified and the ones to be used
            %---------------------------------------------------------------------------------------------------------
            gMarker_Mod      = get(texMarker,'background');
            gMarkerA_Mod     = get(texMarkerA,'background');
            gTissue_Mod      = get(texTissue,'background');
            gTissueA_Mod     = get(texTissueA,'background');
            %------------------------------------------------------------
            g_marker  =           str2num(get(texMarker,'string')); %#ok<*ST2NM>
            g_tissue   =          str2num(get(texTissue, 'string'));
            g_markerA =           str2num(get(texMarkerA,'string'));
            g_tissueA  =          str2num(get(texTissueA, 'string'));
            data.border    =      str2num(get(texBorder,'string'));
            
            
            
            
            %------ if % values change then the abs value needs to be recalculated
            %------ these are the marker values
            %if (floor(1000*g_marker))~= (floor(1000*data.g_marker))
            if gMarker_Mod(1) == 0
                %data.g_markerA =       floor(data.sortFirstSample(floor(data.RCD*g_marker)));   %------ The values are floored for their use in plotting
                data.g_markerA =          (data.sortFirstSample(floor(data.RCD*g_marker)));   %------ No FLOOR for better figures
                data.g_marker =           g_marker;
            elseif gMarkerA_Mod(1) == 0
                data.g_markerA =          g_markerA;
                data.g_marker =           sum(data.sortFirstSample<data.g_markerA)/data.RCD;
            else
                data.g_markerA =          g_markerA;
                data.g_marker =           g_marker;
            end
            %------ these are the tissue values
            %if (floor(1000*g_tissue))~= (floor(1000*data.g_tissue))
            if gTissue_Mod(1) == 0
                %data.g_tissueA =         floor(data.sortFirstSample(floor(data.RCD*g_tissue)));    %------ but perhaps it is not the best option, check later
                data.g_tissueA =          (data.sortFirstSample(floor(data.RCD*g_tissue)));    %------ but perhaps it is not the best option, check later
                data.g_tissue =           g_tissue;
            elseif gTissueA_Mod(1) == 0
                data.g_tissueA =          g_tissueA;
                data.g_tissue =           sum(data.sortFirstSample<data.g_tissueA)/data.RCD;
            else
                data.g_tissue =           g_tissue;
                data.g_tissueA =          g_tissueA;
            end
            set(texTissue, 'background',[1 1 1]);
            set(texTissueA,'background',[1 1 1]);
            set(texMarker, 'background',[1 1 1]);
            set(texMarkerA,'background',[1 1 1]);
            
            
            
            texVoxSizeR =         findobj(gcf,'tag','texVoxSizeR');
            texVoxSizeC =         findobj(gcf,'tag','texVoxSizeC');
            texVoxSizeZ =         findobj(gcf,'tag','texVoxSizeZ');
            
            data.dimVox_R  =           str2num(get(texVoxSizeR,'string'));
            data.dimVox_C  =           str2num(get(texVoxSizeC,'string'));
            data.dimVox_Z  =           str2num(get(texVoxSizeZ,'string'));
            
            set(texVoxSizeR,'background',[1 1 1]);
            set(texVoxSizeC,'background',[1 1 1]);
            set(texVoxSizeZ,'background',[1 1 1]);
            
            %verify that values are within range, relatives [0-tissue-marker-1], abs [0-tis-mark-maxval]
            
            g_tissue_new            = min(max(data.g_tissue,0),data.g_marker);
            g_marker_new            = min(max(data.g_marker,data.g_tissue),1);
            g_tissueA_new           = min(max(data.g_tissueA,0),data.g_markerA);
            g_markerA_new           = min(max(data.g_markerA,data.g_tissueA),data.sortFirstSample(end));
            
            
            if g_tissue_new~=data.g_tissue
                set(texTissue,'string',num2str(g_tissue_new));
                data.g_tissue = g_tissue_new;
            end
            
            if g_marker_new~=data.g_marker
                set(texMarker,'string',num2str(g_marker_new));
                data.g_marker = g_marker_new;
            end
            
            if g_tissueA_new~=data.g_tissueA
                set(texTissueA,'string',num2str(g_tissueA_new));
                data.g_tissueA = g_tissueA_new;
            end
            
            if g_markerA_new~=data.g_markerA
                set(texMarkerA,'string',num2str(g_markerA_new));
                data.g_markerA = g_markerA_new;
            end
  
            
            
            
            %------- recalculate the classes with the new marker values
            %radButtonAdaptHist  = findobj(gcf,'tag','radButtonAdaptHist'); multipleTimeThres = get(radButtonAdaptHist,'value');
%%%%%%%%%%This section was used in case the data was histogram equalised, and several thresholds discard now            
%             if (adaptiveThres == 1)&(multipleTimeThres == 0)
%                 try(data.dataHist);
%                     dataHist = data.dataHist(:,:,:,1);
%                 catch
%                     for counSlice = 1:data.levs
%                         dataTempH = data.dataIn(:,:,counSlice,1);
%                         maxDataTempH = max(dataTempH(:));
%                         dataHist(:,:,counSlice) = maxDataTempH*adapthisteq(dataTempH/maxDataTempH)  ;
%                     end
%                 end
%             elseif (adaptiveThres == 1)&(multipleTimeThres == 1)
%                 try(data.dataHist);
%                     dataHist = data.dataHist;
%                 catch
%                     for counTime = 1:data.timeSamples
%                         for counSlice = 1:data.levs
%                             dataTempH = data.dataIn(:,:,counSlice,counTime);
%                             maxDataTempH = max(dataTempH(:));
%                             dataHist(:,:,counSlice,counTime) = maxDataTempH*adapthisteq(dataTempH/maxDataTempH)  ;
%                         end
%                     end
%                 end
%             elseif (adaptiveThres~= 1)&(multipleTimeThres == 0)
%                 dataHist = data.dataIn(:,:,:,1);
%             elseif (adaptiveThres~= 1)&(multipleTimeThres == 1)
%                 dataHist = data.dataIn(:,:,:,:);
%             end
            dataHist = data.dataIn(:,:,:,1);


            %                [data.classA]  =          getVessels(data.dataIn(:,:,:,1),data.g_markerA,data.g_tissueA,data.border);
            %if (multipleTimeThres == 0)
            %    try (data.timeThresholds);
                    [data.classA]  =    getVessels(dataHist,data.g_markerA,data.g_tissueA,data.border);%[data.classA]  =    getVessels(dataHist,data.timeThresholds(1,1),data.timeThresholds(1,2),data.border);                    
            %    catch
                    
            %    end
            %else
            %    [data.classA]  =          getVessels(dataHist,data.timeThresholds,data.border);
            %end
            %[t1,t2] = contiguity(data.classA(:,:,1),5,4);
            %------- if areas have been manually selected, then remove from classes and append to "other"
            if (sum(sum(data.maskD(:,:,1)))>0)
                %------- once the mask is obtained, remove from the 3 classes:vessels,tissue,border
                data.classA = data.classA.*(1-data.maskD);
            end
            data = getPermeability(data);

            %------ update the objects with the new calculated values
            set(upThresHandle1,'xdata',[data.g_markerA data.g_markerA]);
            set(upThresHandle2,'ydata',[data.g_markerA data.g_markerA]);
            set(lowThresHandle1,'xdata',[data.g_tissueA data.g_tissueA]);
            set(lowThresHandle2,'ydata',[data.g_tissueA data.g_tissueA]);
            set(texMarkerA,'string',num2str(data.g_markerA));
            set(texTissueA,'string',num2str(data.g_tissueA));
            set(texMarker,'string',num2str(data.g_marker,3));
            set(texTissue,'string',num2str(data.g_tissue,3));
            %------ re SURF the classes
            set(classesHandle,  'cdata',(data.classA(1:end,1:end,sliceToPlot)));
            %set(classesHandle,  'cdata',(data.classA(1:end,1:end,sliceToPlot)*4+data.classB(1:end,1:end,sliceToPlot)*3+data.classC(1:end,1:end,sliceToPlot)*2+data.classD(1:end,1:end,sliceToPlot)));
            %------ The measurements in time

            set(lin1,'ydata',data.measurementsInTime(:,1));
            set(lin2,'ydata',data.measurementsInTime(:,2));
            set(lin3,'ydata',data.measurementsInTime(:,3));
            set(lin4,'ydata',data.measurementsInTime(:,4));
            set(lin5,'ydata',data.measurementsInTime(:,5));
            set(linRatio,'ydata',data.measurementsInTime(:,2)./data.measurementsInTime(:,1));
            set(linElements(1),'ydata',data.elementsPerClass(:,5));
            set(linElements(2),'ydata',data.elementsPerClass(:,4));
            set(linElements(3),'ydata',data.elementsPerClass(:,3));
            set(linElements(4),'ydata',data.elementsPerClass(:,2));
            set(linElements(5),'ydata',data.elementsPerClass(:,1));

            set(linPatlak(1),'ydata',data.patlak(:,3),'xdata',data.patlak(:,1));
            set(linPatlak(2),'ydata',data.patlak(:,2),'xdata',data.patlak(:,1));
            set(axTim,'ylimmode','auto');

            set(finVessText,'string', strcat('V(',  num2str(data.timeSamples),') = '  ,num2str(data.measurementsInTime(end,1),2)));
            set(finTissText,'string', strcat('T(',  num2str(data.timeSamples),') = '  ,num2str(data.measurementsInTime(end,2),2)));
            set(finRatioText,'string',strcat('T/V(',num2str(data.timeSamples),') = '  ,num2str(data.measurementsInTime(end,2)./data.measurementsInTime(end,1),2)));

            set(k1Text,'string',strcat('k1 = ',num2str(data.patlak(1,4),3)));
            set(k2Text,'string',strcat('k2 = ',num2str(data.patlak(2,4),3)));
            set(k3Text,'string',strcat('k3 = ',num2str(data.patlak(3,4),3)));
            set(permeabText,'string',strcat('P = ',num2str(data.permeability,3),' [1e7 cm/s]  PS/V = ',num2str(data.PSV,3),' [1e4 1/s]'));

            %set(get(lin1,'parent'),'ylim',[min(min(data.measurementsInTime)) max(max(data.measurementsInTime))] );
            set(fig,'userdata',data);clear
            %------  END of UPDATE ** u_4 **  where RECALCULATION of everything is required --------------------------
        elseif dataIn(3) == '5'
            %---------------------------------------------------------------------------------------------------------
            %------- This is to deselect the ROI
            %---------------------------------------------------------------------------------------------------------
            data.maskD = zeros(data.rows,data.cols,data.levs);
            set(fig,'userdata',data);clear;
            %------  END of UPDATE ** u_5 **   -----------------------------------------------------------------------
        elseif dataIn(3) == '6'
            %---------------------------------------------------------------------------------------------------------
            %------ This is in case a certain rectangular region is the only region of interest
            %---------------------------------------------------------------------------------------------------------
            figure(22)
            imagesc(dataToPlot);
            rect = getrect(figure(22));
            close(figure(22));
            if rect(3)>0
                minRow = max(1,ceil(rect(2)));
                minCol = max(1,ceil(rect(1)));
                maxRow = min(data.rows,floor(rect(2)+rect(4)));
                maxCol = min(data.cols,floor(rect(1)+rect(3)));
            end
            %fig;
            clf;
            PermeabilityAnalysis(data.dataIn(minRow:maxRow,minCol:maxCol,:,:));
            clear;
            %------  END of UPDATE ** u_6 **  where we ZOOM into a ROI -----------------------------------------------
        elseif dataIn(3) == '7'
            %---------------------------------------------------------------------------------------------------------
            %-------- HEre the objects are updated, lines, and values surfs are updated in the next case
            %---------------------------------------------------------------------------------------------------------
            subplot(axisTimeLines);
            switch typeOfPlot
                case 1
                    %----- vessel and tissue ONLY
                    set(lin1,'visible','on');set(lin2,'visible','on');
                    set(lin3,'visible','off');set(lin4,'visible','off');set(lin5,'visible','off');
                    set(linRatio,'visible','off');                    %                set(linPatlak,'visible','off');
                    set(linElements,'visible','off');
                    set(linThresholds,'visible','off');
                    ylimMin = min(min(data.measurementsInTime(:,[1 2])));
                    ylimMax = max(max(data.measurementsInTime(:,[1 2])));
                    xlabel('Time');ylabel('Intensity')
                case 2
                    %----- vessel, tissue and boundaries
                    set(lin1,'visible','on');set(lin2,'visible','on');
                    set(lin3,'visible','on');set(lin4,'visible','on');set(lin5,'visible','off');
                    set(linRatio,'visible','off');%                 set(linPatlak,'visible','off');
                    set(linElements,'visible','off');
                    set(linThresholds,'visible','off');
                    ylimMin = min(min(data.measurementsInTime(:,[1 2 3 4])));
                    ylimMax = max(max(data.measurementsInTime(:,[1 2 3 4])));
                    xlabel('Time');ylabel('Intensity')
                case 3
                    %----- tissue and boundaries plus uncertainty
                    set(lin1,'visible','on');set(lin2,'visible','on');
                    set(lin3,'visible','on');set(lin4,'visible','on');set(lin5,'visible','on');
                    set(linRatio,'visible','off');%                  set(linPatlak,'visible','off');
                    set(linElements,'visible','off');
                    set(linThresholds,'visible','off');
                    ylimMin = min(min(data.measurementsInTime(:,[1 2 3 4 5 ])));
                    ylimMax = max(max(data.measurementsInTime(:,[1 2 3 4 5])));
                    xlabel('Time');ylabel('Intensity')
                case 4
                    %----- boundaries ONLY
                    set(lin1,'visible','off');set(lin2,'visible','on');
                    set(lin3,'visible','on');set(lin4,'visible','on');set(lin5,'visible','off');
                    set(linRatio,'visible','off');%                   set(linPatlak,'visible','off');
                    set(linElements,'visible','off');
                    set(linThresholds,'visible','off');
                    ylimMin = min(min(data.measurementsInTime(:,[2 3 4])));
                    ylimMax = max(max(data.measurementsInTime(:,[2 3 4])));
                    xlabel('Time');ylabel('Intensity')
                case 5
                    %----- vessel, tissue and boundaries +UNCERTAIN
                    set(lin1,'visible','off');set(lin2,'visible','on');
                    set(lin3,'visible','on');set(lin4,'visible','on');set(lin5,'visible','on');
                    set(linRatio,'visible','off');%                    set(linPatlak,'visible','off');
                    set(linElements,'visible','off');
                    set(linThresholds,'visible','off');
                    ylimMin = min(min(data.measurementsInTime(:,[2 3 4 5])));
                    ylimMax = max(max(data.measurementsInTime(:,[2 3 4 5])));
                    xlabel('Time');ylabel('Intensity')
                case 6
                    %------- This is the case of Patlak Plots
                    set(lin1,'visible','off');set(lin2,'visible','off');
                    set(lin3,'visible','off');set(lin4,'visible','off');set(lin5,'visible','off');
                    set(linElements,'visible','off');
                    set(linThresholds,'visible','off');
                    ylimMin = min((min(data.patlak(:,2:3))));
                    ylimMax = max((max(data.patlak(:,2:3))));
                    set(linRatio,'visible','off');%----- Patlak
                    xlimMin = 0.9*(min(data.patlak(:,1)));
                    xlimMax = 1.05*(max(data.patlak(:,1)));
                    set(linPatlak,'visible','on');
                    set(kAText,'visible','on');
                    set(kBText,'visible','off');
                    set(k1Text,'visible','on');
                    set(k2Text,'visible','on');
                    set(k3Text,'visible','on');
                    set(permeabText,'visible','on');
                    set(finVessText,'visible','off');
                    set(finTissText,'visible','off');
                    set(finRatioText,'visible','off');
                    xlabel('Adjusted Time');ylabel('Tissue/Vessel intensity')
                    %-------------------------------------------------------------------------
                case 7
                    set(lin1,'visible','off');set(lin2,'visible','off');set(linPatlak,'visible','off');
                    set(lin3,'visible','off');set(lin4,'visible','off');set(lin5,'visible','off');
                    set(linRatio,'visible','on');
                    set(linElements,'visible','off');
                    set(linThresholds,'visible','off');
                    set(linRatio,'ydata',data.measurementsInTime(:,2)./data.measurementsInTime(:,1));
                    ylimMin = min((min(data.measurementsInTime(:,2)./data.measurementsInTime(:,1)  )));
                    ylimMax = max((max(data.measurementsInTime(:,2)./data.measurementsInTime(:,1)  )));
                    xlabel('Time');ylabel('Intensity RATIO')
                    %----- tissue/vessel intensity
                case 8
                    set(lin1,'visible','on');set(lin2,'visible','on');set(linPatlak,'visible','off');
                    set(lin3,'visible','off');set(lin4,'visible','off');set(lin5,'visible','off');
                    set(linRatio,'visible','off');
                    set(linThresholds,'visible','on');
                    set(linElements,'visible','off');
                    %ylimMin = min((min(data.timeThresholds)));
                    %ylimMax = max((max(data.timeThresholds)));
                                        ylimMin = min(min(data.measurementsInTime(:,[1 2])));
                    ylimMax = max(max(data.measurementsInTime(:,[1 2])));
                    
                    xlabel('Time');ylabel('Intensity')
                    %-----          Intensity thresholds
                case 9
                    set(lin1,'visible','off');set(lin2,'visible','off');set(linPatlak,'visible','off');
                    set(lin3,'visible','off');set(lin4,'visible','off');set(lin5,'visible','off');
                    set(linRatio,'visible','off');
                    set(linThresholds,'visible','off');
                    set(linElements,'visible','on');
                    set(gca,'colororder',[0 0 0.4;0 0.5 0.95;0.4 0.8 0.4 ;1 0.5 0 ;0.5 0 0 ]);
                    ylimMin = min((min(data.elementsPerClass)));
                    ylimMax = max((max(data.elementsPerClass)));
                    xlabel('Time');ylabel('Number of Voxels')
                    %-----         number of voxels in each class
            end
            if typeOfPlot~= 6
                xlimMin = 1;
                xlimMax = data.timeSamples;
                    set(kAText,'visible','off');
                    set(kBText,'visible','on');

                set(linPatlak,'visible','off');
                set(k1Text,'visible','off');
                set(k2Text,'visible','off');
                set(k3Text,'visible','off');
                set(permeabText,'visible','off');
                set(finVessText,'visible','on');
                set(finTissText,'visible','on');
                set(finRatioText,'visible','on');
            end
            set(get(lin1,'parent'),'ylim',[ ylimMin ylimMax],'xlim',[ xlimMin xlimMax] );
            hold off;
            clear;
            %------  END of UPDATE ** u_7 **  where all the OBJECTS (lines, values, ... are updated ----------------------
        elseif dataIn(3) == '8'
            %---------------------------------------------------------------------------------------------------------
            %----- This is the case of the different surf options:
            %     1  Class mask         2 Mask volume Cut       3 Data volume Cut
            %     4  Data time Cut      5 Volume Projection     6 Vessels in 3D
            %     7  Boundary in 3D     8 Slice Montage         9 Volume Montage
            %    10  Histogram
            %-------------------------------------------------------------------
            %---------------------------------------------------------------------------------------------------------
            switch typeOfSurf
                case 1
                    %----- MASK
                    [rowMT,colMT,levMT,timeMT] = size(data.classA);
                    if timeMT == 1
                        set(classesHandle,  'cdata',(data.classA(:,:,sliceToPlot)));
                    else
                        set(classesHandle,  'cdata',(data.classA(:,:,sliceToPlot,timeToPlot)));
                    end
                    set(classesAxis,'xlim',[0.5 data.cols+0.5]);
                    set(XLine2,'xdata',[colToPlot colToPlot]);
                case 2
                    %----- MASK in volume
                    set(classesHandle,  'cdata',(  squeeze(data.classA(:,colToPlot,:))));
                    set(classesAxis,'xlim',[0.5 data.levs+0.5]);
                    set(XLine2,'xdata',[sliceToPlot sliceToPlot]);
                case 3
                    %----- Data in volume
                    set(classesHandle,  'cdata',( squeeze(data.dataIn(:,colToPlot,:,timeToPlot)) ));
                    set(classesAxis,'xlim',[0.5 data.levs+0.5]);
                    set(XLine2,'xdata',[sliceToPlot sliceToPlot]);
                case 4
                    %----- Time Slice
                    set(classesHandle,  'cdata',squeeze(data.dataIn(:,colToPlot,sliceToPlot,:)));
                    set(classesAxis,'xlim',[0.5 data.timeSamples+0.5]);
                    set(XLine2,'xdata',[timeToPlot timeToPlot]);
                    %                    set(classesHandle,'xdata',[1 256])
                case 5
                    %----- Projection
                    set(classesHandle,  'cdata',sum(data.dataIn(:,:,:,timeToPlot),3));
                    set(classesAxis,'xlim',[0.5 data.cols+0.5]);
                    set(XLine2,'xdata',[colToPlot colToPlot]);
                    %case 6
                    %----- Region Leak
                    %set(classesHandle,  'cdata',stretchCMap(data.regLeak(:,:,sliceToPlot)));
                    %set(classesAxis,'xlim',[0.5 data.cols+0.5]);
                    %set(XLine2,'xdata',[colToPlot colToPlot]);
                    %case 7
                    %----- Time Leak
                    %set(classesHandle,  'cdata',stretchCMap(data.timeLeak(:,:,sliceToPlot),1));
                    %%set(classesHandle,  'cdata',(data.timeLeak(:,:,sliceToPlot)));
                    %set(classesAxis,'xlim',[0.5 data.cols+0.5]);
                    %set(XLine2,'xdata',[colToPlot colToPlot]);
                    %case 8
                    %----- Time Leak Normalised
                    %set(classesHandle,  'cdata',stretchCMap(data.timeLeakN(:,:,sliceToPlot),1));
                    %%set(classesHandle,  'cdata',(data.timeLeak(:,:,sliceToPlot)));
                    %set(classesAxis,'xlim',[0.5 data.cols+0.5]);
                    %set(XLine2,'xdata',[colToPlot colToPlot]);
                case 6
                    %----- Vessels in 3D
                    
                    
                    q22=get(gcf,'colormap');
                    figure(fig+1);
                    clf;
                    % if the data is big, subsample
                    if data.RCD>393216
                        subsampleRate =4;
                    else
                        subsampleRate =2;
                        
                    end
                      
                    
                    hold on;
                    [yyy,xxx]                       = meshgrid(1:subsampleRate:data.rows,1:subsampleRate:data.cols);
                    zz                              = ones(size(xxx));
                    for counterZ =1:data.levs;
                        surf(xxx,yyy,(1+data.levs-counterZ)*zz,data.dataIn(1:subsampleRate:end,1:subsampleRate:end,counterZ,timeToPlot),'edgecolor','none')
                    end
                    axis tight
                    view(130,30)
                    rotate3d on;
                    %colorbar;
                    if data.levs>6
                        alpha(0.4);
                    else
                        alpha(0.6);
                    end
                    typeOfSurf2 =0;
                    set(mask_vs_slice,'value',5)
                    colormap(q22)
                    clear q22;
                    zlabel('Z-Position')
                case 7
                    %----- Vessels in 3D
                    figure(fig+1);
                    clf;
                    surfdat(data.classA == 4,'b.'); %surfdat(data.vesVol);
                    typeOfSurf2 =0;
                    set(mask_vs_slice,'value',5)
                case 8
                    %---- =  Inner boundary in 3d
                    figure(fig+1);
                    surfdat(data.classA == 2,'r.');
                    typeOfSurf2 =0;
                    set(mask_vs_slice,'value',5)
                    
                case 9
                    %----- Vessels in @D + t
                    q22=get(gcf,'colormap');
                    figure(fig+1);
                    clf;
                    % if the data is big, subsample
                    if data.RCD>393216
                        subsampleRate =4;
                    else
                        subsampleRate =2;
                    end
                    
                    hold on;
                    [yyy,xxx]                       = meshgrid(1:subsampleRate:data.rows,1:subsampleRate:data.cols);
                    zz                              = ones(size(xxx));
                    for counterZ =1:data.timeSamples;
                        surf(xxx,yyy,(1+data.timeSamples-counterZ)*zz,data.dataIn(1:subsampleRate:end,1:subsampleRate:end,sliceToPlot,counterZ),'edgecolor','none')
                    end
                    axis tight
                    set(gca,'ztick',(1:data.timeSamples),'zticklabel',(data.timeSamples:-1:1))
                    zlabel('Time')

                    view(130,30)
                    rotate3d on;
                    %colorbar;
                    if data.levs>6
                        alpha(0.4);
                    else
                        alpha(0.6);
                    end
                    typeOfSurf2 =0;
                    set(mask_vs_slice,'value',5)
                    colormap(q22)
                    clear q22;

                    
                    
            end
            
            
            switch typeOfSurf2
                case 1
                    %----- a montage of the first volume in a separate figure
                    
                    subplot(plotMontages)
                    hold off
                    dataMontT                   = data.dataIn(:,:,:,timeToPlot);
                    dataMontT(:,end-8:end,:)    = max(dataMontT(:));
                    dataMontT                   = squeeze( dataMontT);
                    [tempRows,tempCols,tempLevs]= size(dataMontT);
                    dataMontT                   = reshape(dataMontT,tempRows,tempCols*tempLevs);
                    imagesc(dataMontT)
                    set(plotMontages,'xtick',-floor(tempCols/2)+linspace(tempCols,tempCols*tempLevs,tempLevs),'xticklabel',linspace(1,tempLevs,tempLevs))
                    
                case 2
                    %----- a montage of the top slice in the handles of the histogram
                    subplot(plotMontages)
                    hold off
                    dataMontT                   = data.dataIn(:,:,sliceToPlot,:);
                    dataMontT(:,end-8:end,:)    = max(dataMontT(:));
                    dataMontT                   = squeeze( dataMontT);
                    [tempRows,tempCols,tempLevs]= size(dataMontT);
                    dataMontT                   = reshape(dataMontT,tempRows,tempCols*tempLevs);
                    imagesc(dataMontT)
                    set(plotMontages,'xtick',-floor(tempCols/2)+linspace(tempCols,tempCols*tempLevs,tempLevs),'xticklabel',linspace(1,tempLevs,tempLevs))
                    
                case 3
                    %----- the histogram of the data
                    
                    subplot(plotMontages)
                    tempHistToPlot          = data.dataIn(:,:,sliceToPlot,timeToPlot);
                    [tempHistY,tempHistX]   = hist(tempHistToPlot(:),100);
                    hold off;
                    bar(data.xHist,data.yHist/sum(data.yHist));
                    hold on;
                    plot(tempHistX,tempHistY/sum(tempHistY),'m')
                    maxYLevel = max(max(data.yHist/sum(data.yHist)),max(tempHistY/sum(tempHistY)));
                    upThresHandle1          = line([data.g_markerA data.g_markerA], [0 maxYLevel],'color',[1 0  0],'tag','upThresHandle1');
                    lowThresHandle1         = line([data.g_tissueA  data.g_tissueA],[0 maxYLevel],'color',[1 0  0],'tag','lowThresHandle1');axis tight
%                    upThresHandle1          = line([data.g_markerA data.g_markerA],[1 data.maxY],'color',[1 0  0],'tag','upThresHandle1');
%                    lowThresHandle1         = line([data.g_tissueA  data.g_tissueA], [1 data.maxY],'color',[1 0  0],'tag','lowThresHandle1');axis tight

            end
            %------  END of UPDATE ** u_8 **  where SURFS (data, masks, ...) are updated-----------------------------
        elseif dataIn(3) == '9'
            %---------------------------------------------------------------------------------------------------------
            %------- REMOVE Time samples
            %------- Read values from a new window
            %---------------------------------------------------------------------------------------------------------
            figure(fig+1);
            fig2 = gcf;
            set(fig2,'MenuBar','none');
            set(fig2,'position',[100 600 600 50])
            set(fig2,'Name',strcat('Time Samples to remove'));
            set(fig2,'NumberTitle','off');
            callBackTimeSamples = ['fig2 = gcf;','fig = fig2-1;','hand3 = get(fig2,''children'');','data = get(fig,''userdata'');',...
                'data.tSeries = str2num(get(hand3,''string''));','close(fig2);',...
                'data.measurementsInTime =  data.measurementsInTime(data.tSeries,:);','data.timeSamples = length(data.tSeries);',...
                'data.elementsPerClass =  data.elementsPerClass(data.tSeries,:);',...
                'data.timeThresholds =  data.timeThresholds(data.tSeries,:);',...
                'data.dataIn = data.dataIn(:,:,:,data.tSeries);','data.patlak = patlak(data.measurementsInTime(:,1:2));',...
                'data.k2                  =  data.patlak(2,4);','data.permeability        =  [1e7*data.k2/60/data.surf_vol];',...
                'set(fig,''userdata'',data);','PermeabilityAnalysis(data);','clear;'];
            texTimeSamples  =  uicontrol(fig2,'Style','edit','string',num2str([1:data.timeSamples]),'position',[10 10 580 25], 'callback',callBackTimeSamples,'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texTimeSamples');
            %texTimeSamples2  =  uicontrol(fig2,'Style','text','string','Please delete time Samples to be removed','position',[10 40 580 25]);
            %------  END of UPDATE ** u_9 **  where TIME SAMPLES are removed -------------------------------------------
        elseif dataIn(3) == 'A'
            %---------------------------------------------------------------------------------------------------------
            %------- REMOVE Data Slices
            %------- Read values from a new window
            %---------------------------------------------------------------------------------------------------------
            figure(fig+1);
            fig2 = gcf;
            set(fig2,'MenuBar','none');
            set(fig2,'position',[100 600 600 50])
            set(fig2,'Name',strcat('Data Slices to remove'));
            set(fig2,'NumberTitle','off');
            callBackSlices = ['fig2 = gcf;','fig = fig2-1;','hand3 = get(fig2,''children'');','data = get(fig,''userdata'');',...
                'data.tSeries = str2num(get(hand3,''string''));',...
                'close(fig2);',...
                'PermeabilityAnalysis(data.dataIn(:,:,data.tSeries,:));','clear;'];
            %                'data.dataIn = data.dataIn(:,:,data.tSeries,:);','[data.rows,data.cols,data.levs,data.timeSamples] = size(data.dataIn);',...
            %                'data.RCD = data.rows*data.cols*data.levs;','set(fig,''userdata'',data);',...
            %                'PermeabilityAnalysis(data);','clear;'];
            texTimeSamples  =  uicontrol(fig2,'Style','edit','string',num2str([1:data.levs]),'position',[10 10 580 25], 'callback',callBackSlices,'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texTimeSamples');
            %texTimeSamples2  =  uicontrol(fig2,'Style','text','string','Please delete time Samples to be removed','position',[10 40 580 25]);
            %------  END of UPDATE ** u_A **  where DATA SLICES are removed -------------------------------------------
        elseif dataIn(3) == 'B'
            %---------------------------------------------------------------------------------------------------------
            %------- Zoom into the data
            %------- Read values from a new window
            %---------------------------------------------------------------------------------------------------------
            figure(fig+1);
            fig2 = gcf;
            set(fig2,'MenuBar','none');
            set(fig2,'position',[100 500 180 200])
            set(fig2,'Name',strcat('Zoom into data'));
            set(fig2,'NumberTitle','off');
            callBackZoomIn = ['fig2 = gcf;','fig = fig2-1;','hand3 = get(fig2,''children'');','data = get(fig,''userdata'');',...
                'newRowI = str2num(get(findobj(hand3,''tag'',''texRowInit''),''string''));',...
                'newRowF = str2num(get(findobj(hand3,''tag'',''texRowFin''),''string''));',...
                'newColI = str2num(get(findobj(hand3,''tag'',''texColInit''),''string''));',...
                'newColF = str2num(get(findobj(hand3,''tag'',''texColFin''),''string''));',...
                'close(fig2);',...
                'data2 = data.dataIn(newRowI:newRowF,newColI:newColF,:,:);',...
                'PermeabilityAnalysis(data2);','clear;'];
            texRowInit   =  uicontrol(fig2,'Style','edit','string','1','position',                [90 140 50 25],'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texRowInit');
            texRowInit2  =  uicontrol(fig2,'Style','text','string','Initial Row','position',      [10 140 70 25],'Handlevisibility','on','BackgroundColor',0.7*[1,1,1],'tag','texRowInit2');
            texRowFin   =  uicontrol(fig2,'Style','edit','string',num2str(data.rows),'position',  [90 110 50 25],'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texRowFin');
            texRowFin2  =  uicontrol(fig2,'Style','text','string','Final Row','position',         [10 110 70 25],'Handlevisibility','on','BackgroundColor',0.7*[1,1,1],'tag','texRowFin2');
            texColInit   =  uicontrol(fig2,'Style','edit','string','1','position',                [90 80 50 25],'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texColInit');
            texColInit2  =  uicontrol(fig2,'Style','text','string','Initial Column','position',   [10 80 70 25],'Handlevisibility','on','BackgroundColor',0.7*[1,1,1],'tag','texColInit2');
            texColFin   =  uicontrol(fig2,'Style','edit','string',num2str(data.cols),'position',  [90 50 50 25],'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texColFin');
            texColFin2  =  uicontrol(fig2,'Style','text','string','Final Column','position',      [10 50 70 25],'Handlevisibility','on','BackgroundColor',0.7*[1,1,1],'tag','texColFin2');


            texZoomCords  =  uicontrol(fig2,'Style','push','string','Zoom','position',[10 10 140 35], 'callback',callBackZoomIn,'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texTimeSamples');
            %texTimeSamples2  =  uicontrol(fig2,'Style','text','string','Please delete time Samples to be removed','position',[10 40 580 25]);
            %------  END of UPDATE ** u_B **  where ZOOM is done through coordinates -------------------------------------------

        elseif dataIn(3) == 'C'
            %---------------------------------------------------------------------------------------------------------
            %------- This is the RE THRESHOLDING  case, recalculate the masks when marker levels are modified FOR THE CURRENT TIME SAMPLE
            %------ Update the values with those from the GUI, the ones with blue background [0 1 1] are the
            %------ last ones to be modified and the ones to be used
            %---------------------------------------------------------------------------------------------------------
            gMarker_Mod             = get(texMarker,'background');
            gMarkerA_Mod            = get(texMarkerA,'background');
            gTissue_Mod             = get(texTissue,'background');
            gTissueA_Mod            = get(texTissueA,'background');
            %------------------------------------------------------------
            g_marker                = str2num(get(texMarker,'string'));
            g_tissue                = str2num(get(texTissue, 'string'));
            g_markerA               = str2num(get(texMarkerA,'string'));
            g_tissueA               = str2num(get(texTissueA, 'string'));
            data.border             = str2num(get(texBorder,'string'));
            %------ if % values change then the abs value needs to be recalculated
            %------ these are the marker values
            %if (floor(1000*g_marker))~= (floor(1000*data.g_marker))
            if gMarker_Mod(1) == 0
                %data.g_markerA     = floor(data.sortFirstSample(floor(data.RCD*g_marker)));   %------ The values are floored for their use in plotting
                data.g_markerA      = (data.sortFirstSample(floor(data.RCD*g_marker)));   %------ No FLOOR for better figures
                data.g_marker       = g_marker;
            elseif gMarkerA_Mod(1) == 0
                data.g_markerA      = g_markerA;
                data.g_marker       = sum(data.sortFirstSample<data.g_markerA)/data.RCD;
            else
                data.g_markerA      = g_markerA;
                data.g_marker       = g_marker;
            end
            %------ these are the tissue values
            %if (floor(1000*g_tissue))~= (floor(1000*data.g_tissue))
            if gTissue_Mod(1) == 0
                %data.g_tissueA     = floor(data.sortFirstSample(floor(data.RCD*g_tissue)));    %------ but perhaps it is not the best option, check later
                data.g_tissueA      = (data.sortFirstSample(floor(data.RCD*g_tissue)));    %------ but perhaps it is not the best option, check later
                data.g_tissue       = g_tissue;
            elseif gTissueA_Mod(1) == 0
                data.g_tissueA      = g_tissueA;
                data.g_tissue       = sum(data.sortFirstSample<data.g_tissueA)/data.RCD;
            else
                data.g_tissue       = g_tissue;
                data.g_tissueA      = g_tissueA;
            end
            set(texTissue, 'background',[1 1 1]);
            set(texTissueA,'background',[1 1 1]);
            set(texMarker, 'background',[1 1 1]);
            set(texMarkerA,'background',[1 1 1]);
            %-------- change the line with the time Thresholds
            data.timeThresholds(timeToPlot,:) = [ data.g_markerA data.g_tissueA];
            if typeOfPlot == 8
            set(linThresholds(1),'ydata',data.timeThresholds(:,1));
            set(linThresholds(2),'ydata',data.timeThresholds(:,2));
            ylimMin                 = min((min(data.timeThresholds)));
            ylimMax                 = max((max(data.timeThresholds)));
            currAxis                = get(linThresholds,'parent');
            set(currAxis{1},'ylim',[ ylimMin ylimMax] );
            end
            %------- recalculate the classes with the new marker values
            if (adaptiveThres == 1)|(typeOfSurf == 13)
                for counSlice = 1:data.levs
                    dataTempH       = data.dataIn(:,:,counSlice,timeToPlot);
                    maxDataTempH    = max(dataTempH(:));
                    dataHist(:,:,counSlice) = maxDataTempH*adapthisteq(dataTempH/maxDataTempH)  ;
                end
                [data.classA]  =          getVessels(dataHist,data.g_markerA,data.g_tissueA,data.border);
            else
                [data.classA]  =          getVessels(data.dataIn(:,:,:,timeToPlot),data.g_markerA,data.g_tissueA,data.border);
            end

            %------- if areas have been manually selected, then remove from classes and append to "other"
            if (sum(sum(data.maskD(:,:,1)))>0)
                %------- once the mask is obtained, remove from the 3 classes:vessels,tissue,border
                data.classA = data.classA.*(1-data.maskD);
            end

            %------ update the objects with the new calculated values
            set(fig,'userdata',data);
            set(upThresHandle1,'xdata',[data.g_markerA data.g_markerA]);
            set(upThresHandle2,'ydata',[data.g_markerA data.g_markerA]);
            set(lowThresHandle1,'xdata',[data.g_tissueA data.g_tissueA]);
            set(lowThresHandle2,'ydata',[data.g_tissueA data.g_tissueA]);
            set(texMarkerA,'string',num2str(data.g_markerA));
            set(texTissueA,'string',num2str(data.g_tissueA));
            set(texMarker,'string',num2str(data.g_marker,3));
            set(texTissue,'string',num2str(data.g_tissue,3));
            %------ re SURF the classes
            set(classesHandle,  'cdata',(data.classA(1:end,1:end,sliceToPlot)));
            %------  END of UPDATE ** u_C **  where RECALCULATION of everything is required --------------------------

        end
    end

end
%-------------End of the updating cases-----------------------------------------------------






%-----------------------------------------------------------------------------------------------------------------------
%---------------------------Extra functions                   ----------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------


%------------------------------------------------------------------
%--------- This is the calculating part of the function ------------
%------------------------------------------------------------------
function   data = calculateVesselData(dataIn,fileName)
%------ it covers TWO cases: (a) the name of the file is received (b) a 4D matrix is received
if ischar(dataIn)
    %------ the received input is a string with the name of the file to process it can be one of two options
    %------              the file with the 4D matrix stored as "B" or
    %------              the file with the processed struct stored as "data"
    %------- keep the name
    load(dataIn);
    if exist('B','var')
        %data.dataIn = B;
        %clear B;
        PermeabilityAnalysis(B);
    else
        data.datName = dataIn;
        PermeabilityAnalysis(data);

    end
end
if ~ischar(dataIn)
    %----------------------------------------------------------------------------------------------------------
    %------ A 4D matrix is received, calculate all necessary values from the input data -----------------------
    %----------------------------------------------------------------------------------------------------------
    data.dataIn                                             = dataIn; 
    clear dataIn;
    if ~exist('fileName','var')
        data.datName                                        = 'data set A';                         %------- No name is received
    else
        data.datName                                        = fileName;
    end

    [data.rows,data.cols,data.levs,data.timeSamples]        = size(data.dataIn);                    %------- Dimensions of data
    data.RCD                                                = data.rows*data.cols*data.levs ;       %------- RCD as total size of elements
    figTemp                                                 = gcf;
    tempData                                                = get(figTemp,'userdata');
    if isstruct(tempData)
        data.g_marker                                       = tempData.g_marker;
        data.g_tissue                                       = tempData.g_tissue;
        data.border                                         = tempData.border;
        clear figTemp tempData
    else

        data.g_marker = 0.95; data.g_tissue = 0.50;  data.border = 5;                               %------- marker values
    end
    data.sortFirstSample                                    = sort(reshape(data.dataIn(:,:,:,1),data.RCD,1));         %------- get the order statistics of the data
    data.g_markerA                                          = min(data.sortFirstSample(end)-1,data.sortFirstSample(floor(data.g_marker*data.RCD)));    %------- only the higher x% are vessel
    data.g_tissueA                                          = data.sortFirstSample(floor(data.g_tissue*data.RCD));    %------- only the lower y% is tissue
    data.colToPlot                                          = floor(data.cols/2);

    %data.dataIn =                 shiftData(data.dataIn,'rigid');                                  %------- shift the data to minimise the movement with time
    tempHistData                                            = data.dataIn(:,:,1,1);
    if ~isa(tempHistData,'double')
        [data.yHist,data.xHist]                                 = hist(double(tempHistData(:)),100);     
    else
        [data.yHist,data.xHist]                                 = hist(tempHistData(:),100); 
    end
    data.maxY                                               = max(data.yHist);                      %------- Histogram of the data for plotting purposes

    %Dimensions of the Voxel
    data.dimVox_R                                           = 1; %[um]  
    data.dimVox_C                                           = 1; %[um]  
    data.dimVox_Z                                           = 1; %[um]  
    
    %------ Segment the data according to the thresholds predefined above----------------------
    [data.classA]                                           =  getVessels(data.dataIn(:,:,:,1),data.g_markerA,data.g_tissueA,data.border);
    %[data.classA,data.classB,data.classC,data.classD]      =  getVessels(data.dataIn(:,:,:,1),data.g_markerA,data.g_tissueA,data.border);
    data.maskD                                              = zeros(data.rows,data.cols,data.levs);
    %------- getPermeability calculates all the below lines in one function
    data                                                    = getPermeability(data);
    % %     %------ With the masks from the segmentation and the data, calculate
    % %     %------             the average value of the class in time
    % %     %------             the patlak adjusted values
    % %     %------             the regional and time leaks
    % %     %------             surface / volume / permeability


    %------------------------------------------------------------------
    %if nargin == 1
        createVesselFigure(data);
    %end
end


%------------------------------------------------------------------
%--------- This is the graphical part of the function ------------
%------------------------------------------------------------------
function   data = createVesselFigure(data)
fig = gcf;

set(fig,'userdata',data);
clf;
scrsz  =  get(0,'ScreenSize');
%Defines window for algorithm and input values
set(fig,'Position',[0.05*scrsz(3) 0.20*scrsz(4) 0.89*scrsz(3) 0.75*scrsz(4)]);   %[scrsz(3)/4 scrsz(4)/4 0.65*scrsz(3) 0.65*scrsz(4)]);
set(fig,'MenuBar','none');
if ~isempty(data)
    set(fig,'Name',strcat('Intensity values in time for: ',data.datName));
else
    set(fig,'Name',strcat('No data has been received '));
end
set(fig,'NumberTitle','off');
set(fig,'Tag','GUIVESSELS');
set(fig,'Toolbar','none','Tag','GUIVESSELS');


%------ File   menu -------------------------------------------------
h_uimenu  =  uimenu(fig,'Label','File');
h_uimenu_read  =  uimenu(h_uimenu,'Label','Read (mat format)','Callback','[sfile,pathname] = uigetfile(''*.mat'',''Please Indicate file'');if sfile~= 0 PermeabilityAnalysis(strcat(pathname,sfile)); end;');
%h_uimenu_read2  =  uimenu(h_uimenu,'Label','Read (biorad format)','Callback','[sfile,pathname] = uigetfile(''*.PIC'',''Please Indicate file'',''multiselect'',''on'');if size(sfile,2)>= 1 data = readPICs(pathname,sfile);PermeabilityAnalysis(data); end;');
h_uimenu_read2  =  uimenu(h_uimenu,'Label','Read (biorad or images)','Callback',' dataIn=readPermeabilityData();PermeabilityAnalysis(dataIn);');
h_uimenu_save  =  uimenu(h_uimenu,'Label','Save','Callback',...
    ['data = get(gcf,''userdata'');[sfile,pathname] = uiputfile(''*.mat'',''Save file as'',strcat(data.datName,''.mat'')); if sfile~= 0 save(strcat(pathname,sfile),''data''); end;']);
h_uimenu_close  =  uimenu(h_uimenu,'Label','Close','Callback',...
    ['close_conf = questdlg(''Do you want to exit?'',''Confirmation'',''Yes'',''No'',''No'');','if strcmp(close_conf,''Yes'');','close;','end;']);

%------ Colour menu -------------------------------------------------
h_uimenu2  =  uimenu(fig,'Label','Colour Map');
h_uimenu2_jet  =      uimenu(h_uimenu2,'Label','jet',     'Callback','colormap(jet);');
h_uimenu2_jet2      =   uimenu(h_uimenu2,'Label','jet (low)', 'Callback','qq(:,3) =[linspace(0.5,1,8)'';ones(16,1);linspace(1,0,15)'';zeros(25,1)] ;qq(:,1) = [zeros(24,1);linspace(0,1,15)'';ones(16,1);linspace(1,0.5,9)''];qq(:,2) = [zeros(8,1);linspace(0,1,15)'';ones(16,1);linspace(1,0,15)'';zeros(10,1)]; colormap(qq([1:3:40 41:end],:));clear qq;');
h_uimenu2_jet3      =   uimenu(h_uimenu2,'Label','jet (high)','Callback','qq(:,3) =[linspace(0.5,1,8)'';ones(16,1);linspace(1,0,15)'';zeros(25,1)] ;qq(:,1) = [zeros(24,1);linspace(0,1,15)'';ones(16,1);linspace(1,0.5,9)''];qq(:,2) = [zeros(8,1);linspace(0,1,15)'';ones(16,1);linspace(1,0,15)'';zeros(10,1)]; colormap(qq([1:20 21:2:end],:));clear qq;');


h_uimenu2_hot  =      uimenu(h_uimenu2,'Label','hot',     'Callback','colormap(hot);');
h_uimenu2_hot2  =     uimenu(h_uimenu2,'Label','hot^2',   'Callback','colormap(hot.^2);');
h_uimenu2_hot12  =     uimenu(h_uimenu2,'Label','hot^(1/2)','Callback','colormap(hot.^(0.5));');

h_uimenu2_gray  =     uimenu(h_uimenu2,'Label','gray',    'Callback','colormap(gray);');
h_uimenu2_gray2  =    uimenu(h_uimenu2,'Label','gray^(1/2)','Callback','colormap(gray.^(0.5));');
h_uimenu2_gray12  =    uimenu(h_uimenu2,'Label','gray^2',  'Callback','colormap(gray.^2);');

%------ Algorithm menu -------------------------------------------------
% h_uimenu3  =  uimenu(fig,'Label','Algorithms');
% h_uimenu3_reduce  =      uimenu(h_uimenu3,'Label','** Reduce size (Quad/Oct Tree)',               'Callback','fig = gcf;data = get(fig,''userdata'');data2 = reduceu(data.dataIn);PermeabilityAnalysis(data2);clear;');
% h_uimenu3_rotate  =      uimenu(h_uimenu3,'Label','Transpose data (x<--->y)',                     'Callback','fig = gcf;data = get(fig,''userdata'');data.classA = permute(data.classA,[2 1 3 4]);data.maskD = permute(data.maskD,[2 1 3 4]);data.dataIn = permute(data.dataIn,[2 1 3 4]);temRow = data.rows;data.rows = data.cols;data.cols = temRow;PermeabilityAnalysis(data);clear;');
% h_uimenu3_register1  =   uimenu(h_uimenu3,'Label','** Registration of data (Rigid)',              'Callback','fig = gcf;data = get(fig,''userdata'');data2 = shiftData(data.dataIn,''rigid'');PermeabilityAnalysis(data2);clear;');
% h_uimenu3_register2  =   uimenu(h_uimenu3,'Label','** Registration of data (Regional Rigid)',     'Callback','fig = gcf;data = get(fig,''userdata'');data2 = shiftData(data.dataIn,''nonrigid'');PermeabilityAnalysis(data2);clear;');
% h_uimenu3_zoom  =        uimenu(h_uimenu3,'Label','** Zoom into a region (graph) ',               'Callback','PermeabilityAnalysis(''u_6'');');
% h_uimenu3_zoom2  =       uimenu(h_uimenu3,'Label','** Zoom into a region (coords)',               'Callback','PermeabilityAnalysis(''u_B'');');
% h_uimenu3_equalise  =    uimenu(h_uimenu3,'Label','** Equalise intensity in the slides',          'Callback','fig = gcf;data = get(fig,''userdata'');avIntensity = (sum(sum(data.dataIn(:,:,:,1)))/data.rows/data.cols);avIntensityDiff = avIntensity-avIntensity(1); q = repmat(avIntensityDiff,[data.rows data.cols 1 data.timeSamples]);   data.dataIn = data.dataIn-q;figure(fig+1);plot(squeeze(avIntensity),''b-o'');figure(fig);PermeabilityAnalysis(data);clear;');
% h_uimenu3_warn  =        uimenu(h_uimenu3,'Label','** Warning! these process ** cannot be undone ');


if ~isempty(data)
    %%
    BackFrame           =  uicontrol(fig,'Style','frame','position',[5  5 580   190],    'visible','on','tag','BackFrame');
  %%
    % ------ Push Button to End programme
    pbend =  uicontrol(fig,'Style','push','String','End',                                 'position',[10 15 150 30],...
        'CallBack','data = get(gcf,''userdata'');close;return;','Handlevisibility','on','tag','pbend','value',0);
    %------ Push Button to Recalculate and Replot everything
    pbgraf =  uicontrol(fig,'Style','push','String','Recalculate',                        'position',[10 45 150 30],...
        'CallBack','PermeabilityAnalysis(''u_4'')','Handlevisibility','on','tag','pbgraf','value',0);
    %------ Push Button to  to remove a  mask
    %pbmask =  uicontrol(fig,'Style','push','String','Erase ROI',                          'position',[170 35 110 30],...
    %    'CallBack','PermeabilityAnalysis(''u_3'')','Handlevisibility','on','tag','pbmask','value',0);
    %------ Push Button to  Zoom int a Region of Interest
    %pbzoom =  uicontrol(fig,'Style','push','String','Zoom ROI',                           'position',[285 35 110 30],...
    %    'CallBack','PermeabilityAnalysis(''u_6'')','Handlevisibility','on','tag','pbzoom','value',0);
    %------ button do a thresholding of the current time sample
    %pbThreshold  =  uicontrol(fig,'Style','push','String','Threshold',                   'position',[170  72 110 30],...
    %    'CallBack','PermeabilityAnalysis(''u_C'')','Handlevisibility','on','tag','pbThres','value',0);
    %------ Push Button to remove deselect the Region of interest
    %pbincludeMask =  uicontrol(fig,'Style','push','String','De-select ROI',               'position',[400 35 110 30],...
    %    'CallBack','PermeabilityAnalysis(''u_5'')','Handlevisibility','on','tag','pbinclude','value',0,'fontsize',10);
    %------ button to remove time samples
    pbTSamps =  uicontrol(fig,'Style','push','String','Remove samples',                   'position',[170 45 150 30],...
        'CallBack','PermeabilityAnalysis(''u_9'')','Handlevisibility','on','tag','pbTSamps','value',0);

    %------ button to remove SLICES
    pbSlices =  uicontrol(fig,'Style','push','String','Remove Slices',                   'position',[170 15 150 30],...
        'CallBack','PermeabilityAnalysis(''u_A'')','Handlevisibility','on','tag','pbSlices','value',0);
%
    %------  VAlues for the thresholds editable and in %
          uicontrol(fig,'Style','edit','string',num2str(data.g_marker,3),       'position',[10 165 45 22], ...
        'callback',['texMarker =  findobj(gcf,''tag'',''texMarker'');texMarkerA =  findobj(gcf,''tag'',''texMarkerA'');set(texMarker,''background'',[0 1 1]);set(texMarkerA,''background'',[1 1 1]);clear texM*;'],...
        'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texMarker');
          uicontrol(fig,'Style','edit','string',num2str(data.g_tissue,3),       'position',[10 140 45 22], ...
        'callback',('texTissue =  findobj(gcf,''tag'',''texTissue'');texTissueA =  findobj(gcf,''tag'',''texTissueA'');set(texTissue,''background'',[0 1 1]);set(texTissueA,''background'',[1 1 1]);clear texT*;'),...
        'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texTissue');
           uicontrol(fig,'Style','edit','string',num2str(data.colToPlot),        'position',[10 115 45 22],...
        'callback','PermeabilityAnalysis(''u_1'')','Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texXLine');
          uicontrol(fig,'Style','edit','string',num2str(data.border),           'position',[10 90 45 22],...
        'callback','texBorder =  findobj(gcf,''tag'',''texBorder'');set(texBorder,''background'',[0 1 1]);clear texB*;','Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texBorder');

    %------  VAlues for the thresholds, non editable and in pix intensity
    texMarkerA              =  uicontrol(fig,'Style','edit','string',num2str(data.g_markerA),                               'position',[110 165 45 22], ...
        'callback',['texMarker =  findobj(gcf,''tag'',''texMarker'');texMarkerA =  findobj(gcf,''tag'',''texMarkerA'');set(texMarker,''background'',[1 1 1]);set(texMarkerA,''background'',[0 1 1]);clear texM*;'], ...
        'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texMarkerA');
    texTissueA              =  uicontrol(fig,'Style','edit','string',num2str(data.g_tissueA),                               'position',[110 140 45 22], ...
        'callback',['texTissue =  findobj(gcf,''tag'',''texTissue'');texTissueA =  findobj(gcf,''tag'',''texTissueA'');set(texTissue,''background'',[1 1 1]);set(texTissueA,''background'',[0 1 1]);clear texT*;'], ...
        'Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texTissueA');
    %------  Text for the names
    nomMarker               =  uicontrol(fig,'Style','text','string','Marker',                                              'position',[55 140 45 22], 'Handlevisibility','on','tag','nomMarker');
    nomTissue               =  uicontrol(fig,'Style','text','string','Tissue',                                              'position',[55 165 45 22], 'Handlevisibility','on','tag','nomTissue');
    nomBorder               =  uicontrol(fig,'Style','text','string','Border',                                              'position',[55 90  45 22], 'Handlevisibility','on','tag','nomBorder');
    nomXLine                =  uicontrol(fig,'Style','text','string','X Line',                                              'position',[55 115 45 22], 'Handlevisibility','on','tag','nomXLine');

    %-------------------------Change the slice that is being viewed --------------------
    %------  Text for the names
    nomSlice                =  uicontrol(fig,'Style','text','string',strcat('slice: (1:',num2str(data.levs),')'),           'position',[225 140 85 22], 'Handlevisibility','on','tag','nomSlice');
    nomTime                 =  uicontrol(fig,'Style','text','string',strcat('Time:  (1:',num2str(data.timeSamples),')'),    'position',[225 165 85 22], 'Handlevisibility','on','tag','nomTime');
    nomRows                 =  uicontrol(fig,'Style','text','string',strcat('Rows:     ',num2str(data.rows)),               'position',[225 115 85 22], 'Handlevisibility','on','tag','nomSlice');
    nomCols                 =  uicontrol(fig,'Style','text','string',strcat('Columns: ',num2str(data.cols)),                'position',[225 90  85 22], 'Handlevisibility','on','tag','nomSlice');
    texTime                 =  uicontrol(fig,'Style','edit','string','1',                                                   'position',[175 165 45 22], 'callback','PermeabilityAnalysis(''u_2'')','Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texTime');
    texSlice                =  uicontrol(fig,'Style','edit','string','1',                                                   'position',[175 140 45 22], 'callback','PermeabilityAnalysis(''u_2'')','Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texSlice');

    %%
    %-------------------------Change dimensions of the voxel if necessary --------------------
    %------  Text for the names
    %text(355,85,'\mu','interpreter','tex')

    nomVoxSize              =  uicontrol(fig,'Style','text','string','Voxel dimensions',                                    'position',[505 22 70 45], 'Handlevisibility','on','tag','nomVoxSize');
    nomVoxSizeR             =  uicontrol(fig,'Style','text','string','Row     : ',                                          'position',[330 60 70 22], 'Handlevisibility','on','tag','nomVoxSizeR');
    nomVoxSizeC             =  uicontrol(fig,'Style','text','string','Column  : ',                                          'position',[330 40 70 22], 'Handlevisibility','on','tag','nomVoxSizeC');
    nomVoxSizeZ             =  uicontrol(fig,'Style','text','string','Z Level : ',                                          'position',[330 20 70 22], 'Handlevisibility','on','tag','nomVoxSizeZ');

    texVoxSizeR             =  uicontrol(fig,'Style','edit','string',num2str(data.dimVox_R),                                'position',[405 60 50 22], 'callback','texVoxSizeR =  findobj(gcf,''tag'',''texVoxSizeR'');set(texVoxSizeR,''background'',[0 1 1]);clear texV*;','Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texVoxSizeR');
    texVoxSizeC             =  uicontrol(fig,'Style','edit','string',num2str(data.dimVox_C),                                'position',[405 40 50 22], 'callback','texVoxSizeC =  findobj(gcf,''tag'',''texVoxSizeC'');set(texVoxSizeC,''background'',[0 1 1]);clear texV*;','Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texVoxSizeC');
    texVoxSizeZ             =  uicontrol(fig,'Style','edit','string',num2str(data.dimVox_Z),                                'position',[405 20 50 22], 'callback','texVoxSizeZ =  findobj(gcf,''tag'',''texVoxSizeZ'');set(texVoxSizeZ,''background'',[0 1 1]);clear texV*;','Handlevisibility','on','BackgroundColor',[1,1,1],'tag','texVoxSizeZ');

    nomVoxSizeR2            =  uicontrol(fig,'Style','text','string','[m','fontname','symbol',                               'position',[460 60 12 22], 'Handlevisibility','on','tag','nomVoxSizeR2');
    nomVoxSizeR3            =  uicontrol(fig,'Style','text','string','m]',                                                   'position',[472 60 15 22], 'Handlevisibility','on','tag','nomVoxSizeR3');
    nomVoxSizeC2            =  uicontrol(fig,'Style','text','string','[m','fontname','symbol',                               'position',[460 40 12 22], 'Handlevisibility','on','tag','nomVoxSizeC2');
    nomVoxSizeC3            =  uicontrol(fig,'Style','text','string','m]',                                                   'position',[472 40 15 22], 'Handlevisibility','on','tag','nomVoxSizeC3');
    nomVoxSizeZ2            =  uicontrol(fig,'Style','text','string','[m','fontname','symbol',                               'position',[460 20 12 22], 'Handlevisibility','on','tag','nomVoxSizeZ2');
    nomVoxSizeZ3            =  uicontrol(fig,'Style','text','string','m]',                                                   'position',[472 20 15 22], 'Handlevisibility','on','tag','nomVoxSizeZ3');
 
%%    
    
    t_vs_pat                =  uicontrol(fig,'Style','popup','string','Intensity Vessel,tissue|Intensity Vessel,tissue,boundaries|Intensity Vessel,tissue,boundaries,Uncertain|Intensity Tissue,Boundaries|Intensity Tissue,boundaries,Uncertain|Patlak|Intensity Tissue/Vessel|Vessel, Tissue + Thresholds','value',1,'units','normalized', 'position',[0.6 0.01 0.28 0.05],  'callback','PermeabilityAnalysis(''u_7'')','Handlevisibility','on','tag','t_vs_pat');
%    t_vs_pat                =  uicontrol(fig,'Style','popup','string','Intensity Vessel,tissue|Intensity Vessel,tissue,boundaries|Intensity Vessel,tissue,boundaries,Uncertain|Intensity Tissue,Boundaries|Intensity Tissue,boundaries,Uncertain|Patlak|Intensity Tissue/Vessel|Vessel, Tissue + Thresholds','value',1, 'position',[550 10 260 25],  'callback','PermeabilityAnalysis(''u_7'')','Handlevisibility','on','tag','t_vs_pat','units','normalized');
    %mask_vs_slice  =  uicontrol(fig,'Style','popup','string','Class mask (row,col)|Mask volume Cut (row,slice)|Data volume Cut (row,slice)|Data time Cut (row,t)|Volume Projection (row,col)|Region Leak (row,col)|Time Leak (row,col)|Time Leak Norm(row,col)|Vessels (3D)|Boundary (3D)|1 Vol Montage(x,y,z)|1 Slice Montage(x,y,t)|Histogram|Data Slice Linear Montage|Mask Slice Linear Montage','value',1, 'position',[750 210 405 25],  'callback','PermeabilityAnalysis(''u_8'')','Handlevisibility','on','tag','mask_vs_slice');
%%
%    mask_vs_slice           =  uicontrol(fig,'Style','popup','string','Class mask [rows,columns]|Mask volume Cut [rows,z position]|Data volume Cut [rows, z position]|Data time Cut [rows,time]|Average Intensity Projection [rows,columns]|Z-stack Intensity [rows,columns,z position]|Vessels (3D) [rows,columns,z position]|Boundary (3D) [rows,columns,z position]','value',1, 'position',[750 210 405 25],  'callback','PermeabilityAnalysis(''u_8'')','Handlevisibility','on','tag','mask_vs_slice','units','normalized');
    mask_vs_slice           =  uicontrol(fig,'Style','popup','string','Class mask [rows,columns]|Mask volume Cut [rows,z position]|Data volume Cut [rows, z position]|Data time Cut [rows,time]|Average Intensity Projection [rows,columns]|Z-stack Intensity [rows,columns,z position]|Vessels (3D) [rows,columns,z position]|Boundary (3D) [rows,columns,z position]|3D Time Projection [rows,columns,time]','value',1,'units','normalized', 'position',[0.6 0.3  0.28 0.05],  'callback','PermeabilityAnalysis(''u_8'')','Handlevisibility','on','tag','mask_vs_slice');

%%    
%    hist_vs_montage         =  uicontrol(fig,'Style','popup','string','Montage (x,y,z) of one time stack [rows,z position]|Montage (x,y,t) of one slice along time [rows,time] |Histogram  [relative occurrence, voxel intensity]','value',3,'units','normalized', 'position',[750 505 405 25],  'callback','PermeabilityAnalysis(''u_8'')','Handlevisibility','on','tag','hist_vs_montage');
    hist_vs_montage         =  uicontrol(fig,'Style','popup','string','Montage (x,y,z) of one time stack [rows,z position]|Montage (x,y,t) of one slice along time [rows,time] |Histogram  [relative occurrence, voxel intensity]','value',3,'units','normalized', 'position', [0.6 0.756  0.28 0.05],  'callback','PermeabilityAnalysis(''u_8'')','Handlevisibility','on','tag','hist_vs_montage');

%%    
    
   
    %------- this is the histogram and the threshold
    plotMontages            = subplot('Position', [0.53      0.83    0.45    0.15]);
    set(plotMontages,'tag','plotMontages');
    hold off;
    bar(data.xHist,data.yHist/sum(data.yHist));
    hold on;
    
    
    maxYLevel = max(data.yHist/sum(data.yHist));
    upThresHandle1          = line([data.g_markerA data.g_markerA], [0 maxYLevel],'color',[1 0  0],'tag','upThresHandle1');
    lowThresHandle1         = line([data.g_tissueA  data.g_tissueA],[0 maxYLevel],'color',[1 0  0],'tag','lowThresHandle1');axis tight
    
%     upThresHandle1          = line([data.g_markerA data.g_markerA],[1 data.maxY],'color',[1 0  0],'tag','upThresHandle1');
%     lowThresHandle1         = line([data.g_tissueA  data.g_tissueA], [1 data.maxY],'color',[1 0  0],'tag','lowThresHandle1');axis tight
%     
    
    
%%    

    
    
    %------ this is the line plot + thresholds
    %plotHandle             = subplot('Position', [0.56    0.82    0.4    0.16]);
    plotHandle              = subplot('Position', [0.03      0.83    0.4    0.15]);
    set(plotHandle,'tag','plotHandle');
    hold off;
    XLineP                  = plot(1:data.rows,data.dataIn(:,min(data.colToPlot,data.cols),1,1),'b','tag','XLineP');hold on
    XLineP2                 = plot(1:data.rows,data.dataIn(:,min(data.colToPlot,data.cols),1,1),'m','tag','XLineP2','visible','off');hold on
    upThresHandle2          = line([1 data.rows],[data.g_markerA data.g_markerA],'color',[ 1 0 0],'tag','upThresHandle2');
    lowThresHandle2         = line([1 data.rows],[data.g_tissueA  data.g_tissueA],'color',[1 0  0],'tag','lowThresHandle2');axis tight
%
%     dataMontT                   = data.dataIn(1:1:end,1:1:end,1,:);
%     dataMontT(:,end-8:end,:)    = max(dataMontT(:));
%     dataMontT                   = squeeze( dataMontT);
%     [tempRows,tempCols,tempLevs]= size(dataMontT);
%     dataMontT                   = reshape(dataMontT,tempRows,tempCols*tempLevs);
%     imagesc(dataMontT)
%     set(plotHandle,'xtick',-floor(tempCols/2)+linspace(tempCols,tempCols*tempLevs,tempLevs),'xticklabel',linspace(1,tempLevs,tempLevs))

%    
    
    %------ this is the montage plot 
%     plotMontages         = subplot('Position', [0.5    0.83    0.48    0.15]);
%     dataTemp             = [squeeze(data.dataIn(1:2:end,1:2:end,1,:)) 255*ones(data.rows/2,1,data.timeSamples) zeros(data.rows/2,1,data.timeSamples)];
%     dataTemp             = reshape(dataTemp,[data.rows/2 (data.cols/2+2)*data.timeSamples]);
%     imagesc(dataTemp)
    %------ Surf of one slice
                            subplot('Position', [0.03   0.40    0.45    0.35]);
    dataHandle           =  imagesc(data.dataIn(1:end,1:end,1,1)); set(dataHandle,'tag','dataHandle');
    colorbar;
    XLine                =  line([data.colToPlot data.colToPlot],[1 data.rows],'color',[1 1 1],'linestyle','--','tag','XLine');

    %------ the classes
    classesAxis          =  subplot('Position', [0.53    0.40    0.45    0.35]);
    classesHandle        =  imagesc(data.classA(1:end,1:end,1));
    set(classesHandle,'tag','classesHandle');
    set(classesAxis,'tag','classesAxis');
    XLine2               =  line([data.colToPlot data.colToPlot],[1 data.rows],'color',[0 0 0],'linestyle','--','tag','XLine2');
    
    %------ The measurements in time  0.48  0.12 0.48 0.1957
    axisTimeLines        = subplot('Position', [0.53    0.13    0.45    0.16]);
    set(axisTimeLines,'ylimmode','auto','tag','axisTimeLines');
    hold on;
%%
%     %----------------------------------- This is for single/multiple time thresholding
%     radButtonSingle      = uicontrol(fig,'Style','radio','String','Single Time Thresholding',  'position',[55 205 180 25],'value',1,'tag','radButtonSingle' ,'callback','radButtonMulti  = findobj(gcf,''tag'',''radButtonMulti''); set(radButtonMulti,''value'',0);');
%     radButtonMulti       = uicontrol(fig,'Style','radio','String','Multiple Time Thresholding','position',[55 180 180 25],'value',0,'tag','radButtonMulti','callback','radButtonSingle = findobj(gcf,''tag'',''radButtonSingle'');set(radButtonSingle,''value'',0);t_vs_pat = findobj(gcf,''tag'',''t_vs_pat'');set(t_vs_pat,''value'',8);PermeabilityAnalysis(''u_7'');');
% 
%     %----------------------------------- This is for thresholding with original/AdaptiveHistEq data
%     radButtonOriginal    = uicontrol(fig,'Style','radio','String','Thresholding original data',  'position',[255 205 180 25],'value',1,'tag','radButtonOriginal' ,...
%         'callback','radButtonAdaptHist = findobj(gcf,''tag'',''radButtonAdaptHist'');set(radButtonAdaptHist,''value'',0); mask_vs_slice = findobj(gcf,''tag'',''mask_vs_slice'');set(mask_vs_slice,''value'',1);');
%     radButtonAdaptHist   = uicontrol(fig,'Style','radio','String','Thresholding adaptive data',  'position',[255 180 180 25],'value',0,'tag','radButtonAdaptHist',...
%         'callback','radButtonOriginal  = findobj(gcf,''tag'',''radButtonOriginal''); set(radButtonOriginal,''value'',0);  mask_vs_slice = findobj(gcf,''tag'',''mask_vs_slice'');set(mask_vs_slice,''value'',13);PermeabilityAnalysis(''u_8'');');
    
    %dataPSV = data.PSV;dataElementsPerClass = data.elementsPerClass ;dataTimeThresholds = data.timeThresholds ;
    
    
    %------- Verify that all the fields have been generated
    if ~isfield(data,'PSV')       
        data.PSV                                 =  1e4*data.k2/60;
    end
        
    if ~isfield(data,'elementsPerClass')
        [tempVals,data.elementsPerClass]         =  calculateTimeValues(data.dataIn,data.classA); 
        clear tempVals;
        set(fig,'userdata',data);
    end
        
    if ~isfield(data,'timeThresholds')  
        data.timeThresholds                      =  repmat([data.g_markerA data.g_tissueA],[data.timeSamples 1]);
        set(fig,'userdata',data);
    end
    
%     if ~isfield(data,'jet4')       
%         data.jet4                                =  jet;
%         set(fig,'userdata',data);
%     end


    %------- these are the colour codes to be used for the lines so they match the MASK in JET Colour map
    colour5         =   [0 0 0.4 ];      %   dark blue       colour1  =   [0.6 0.3 0.3];   %---- brown
    colour4         =   [0 0.5 0.95 ];   %   blue            colour2  =   [0.9 0.9 0.0] ;  %---- yellow
    colour3         =   [0.4 0.8 0.4 ];  %   green           colour3  =   [0. 0.8 0.8]  ;  %---- cyan
    colour2         =   [1 0.5 0 ];      %   orange          colour4  =   [0. 0. 0.5] ;    %---- blue
    colour1         =   [0.5 0 0 ];      %   brown           colour5  =    ;
    colourOrder     = [colour1;colour2;colour3;colour4;colour5];
    
    %------- Plot the lines of the intensity or the Patlak plot They are all printed at the same time
    %------- and then they are selected on/off according to the plot selected by the pull down menu
    lin1            = plot(1:data.timeSamples,data.measurementsInTime(:,1),'color',colour1,'marker','.'   ,'tag','lin1');        %------ vessels
    lin2            = plot(1:data.timeSamples,data.measurementsInTime(:,2),'color',colour2,'marker','.'    ,'tag','lin2');       %------ tissue
    lin3            = plot(1:data.timeSamples,data.measurementsInTime(:,3),'color',colour3,'marker','.'   ,'tag','lin3');        %------ border1
    lin4            = plot(1:data.timeSamples,data.measurementsInTime(:,4),'color',colour4,'marker','.'    ,'tag','lin4');       %------ border2
    lin5            = plot(1:data.timeSamples,data.measurementsInTime(:,5),'color',colour5,'marker','.'    ,'tag','lin5');       %------ else
    linPatlak       = plot(data.patlak(:,1),data.patlak(:,2),'r*',data.patlak(:,1),data.patlak(:,3),'b--.','visible','off','tag','linPatlak');
    linRatio        = plot(1:data.timeSamples,data.measurementsInTime(:,2)./data.measurementsInTime(:,1),'color',colour1,'marker','.'    ,'tag','linRatio','visible','off');     %------ RAtio of Tissue/Vessel
    linElements     = plot(repmat([1:data.timeSamples]',[1 5]),data.elementsPerClass, 'marker','.'    ,'tag','linElements','visible','off');     %------);
    linThresholds   = plot(1:data.timeSamples,data.timeThresholds(:,1),'r-x',1:data.timeSamples,data.timeThresholds(:,2),'b-o','tag','linThresholds','visible','off');     %------);
    set(linElements(1),'color',colour1);set(linElements(2),'color',colour2);set(linElements(3),'color',colour3);
    set(linElements(4),'color',colour4);set(linElements(5),'color',colour5);
    
    %------- These are the numerical values that correspond to the plots, last values of the vessel/tissue
    %------- or the patlak calculations
    kAText           =  uicontrol(fig,'Style','text','position',[505  140 70   45],   'string','Patlak Parameters', 'visible','off','tag','kAText');
    k1Text           =  uicontrol(fig,'Style','text','position',[330  164 155  22],  'string',strcat('k1 = '    ,num2str(data.patlak(1,4),3)), 'visible','off','tag','k1Text');
    k2Text           =  uicontrol(fig,'Style','text','position',[330  142 155  22],  'string',strcat('k2 = '    ,num2str(data.patlak(2,4),3)), 'visible','off','tag','k2Text');
    k3Text           =  uicontrol(fig,'Style','text','position',[330  120 155  22],  'string',strcat('k3 = '    ,num2str(data.patlak(3,4),3)), 'visible','off','tag','k3Text');
    permeabText      =  uicontrol(fig,'Style','text','position',[330  90  155  30],  'string',strcat('P = '    ,num2str(data.permeability,4),' [1e7 cm/s]  PS/V = ',num2str(data.PSV,3),' [1e4 1/s]'),'visible','off','tag','permeabText');
    kBText           =  uicontrol(fig,'Style','text','position',[505  140 70   45],   'string','Intensity Parameters', 'visible','on','tag','kBText');

    finVessText      =  uicontrol(fig,'Style','text','position',[330  164 155  22],  'string',strcat('V(',num2str(data.timeSamples),') = '  ,num2str(data.measurementsInTime(end,1),4)), 'visible','on','tag','finVessText');
    finTissText      =  uicontrol(fig,'Style','text','position',[330  142 155  22],  'string',strcat('T(',num2str(data.timeSamples),') = '  ,num2str(data.measurementsInTime(end,2),4)), 'visible','on','tag','finTissText');
    finRatioText     =  uicontrol(fig,'Style','text','position',[330  120 155  22],  'string',strcat('T/V(',num2str(data.timeSamples),') = '  ,num2str(data.measurementsInTime(end,2)./data.measurementsInTime(end,1),4)), 'visible','on','tag','finRatioText');
    axis ([1 data.timeSamples min(min(data.measurementsInTime(:,1:5))) max(max(data.measurementsInTime(:,1:5))) ]);grid on;
    xlabel('Time');ylabel('Intensity')
    hold off;
end
%------------------------------------------------------------------
%--------- End of the graphical output of the function ------------
%------------------------------------------------------------------










