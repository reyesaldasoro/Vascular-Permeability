function [dataIn]=readPermeabilityData(dataInName)

%function [data]=readPermeabilityData(dataInName)
%-------------------------------------------------------------------------------------
%------- VARARGIN   :   {1} data  =  image to be analysed it can be an image, mat file or file name
%-------                    it can also be a folder with more folders of a folder with tiffs
%------- ARGOUT     :   handles will keep the number of time frames and
%-------                the data is saved into a new folder with the extension _MAT_OR where the original
%-------                data in matlab format is stored
%-------------------------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro                       ----------
%------             The University of Sheffield                             ----------
%------             http://carlos-reyes.staff.shef.ac.uk                    ----------
%------  20 Jan 2009                                       ---------------------------
%-------------------------------------------------------------------------------------


%% Parse input

switch nargin
    case 0
        %----- no data received, Open question dialog and pass to next section to analyse
        %
        button                                  = questdlg('No data or figure received','Select Input','Single File (Matlab)',sprintf('Files in Folder (mat,tiff,bmp,PIC,...)'),'Cancel','Cancel');
        
        %
        if strcmp(button(1),'S')
            %------ a single file, can ONLY be a 4D matlab matrix
            [sfile,pathname]                    =  uigetfile('*.*','Please Select data to load (*.mat)');
            if sfile~=  0
                dataInName                      =  strcat(pathname,'/',sfile);
                clear sfile pathname;
                dataFromFile                    = load(dataInName);
                if isfield(dataFromFile,'dataIn')
                    dataIn                      = dataFromFile.dataIn;
                else
                    namesF                      = fieldnames(dataFromFile);
                    dataIn                      = getfield(dataFromFile,namesF{1});
                end
                %handles.numFrames               = size(dataIn,4);
            else
                disp('File not found');
                dataIn=[];
                return;
            end
        elseif strcmp(button(1),'M')
            %----- a folder, can be a folder with matlab files, or a folder with folders and tiff images
            [pathname]                          =  uigetdir('*.*','Please select folder where the images/data are');
            if pathname~=  0
                % pass the pathname to same function to process
                [dataIn]                = readPermeabilityData(pathname);
            else
                disp('Folder not found');
                dataIn=[];
                return;
            end
        else
            %disp('Cancel');
            dataIn=[];
            return;
        end
    case 1
        %----- one argument received, it can be a char of a file name or a folder
        if isa(dataInName,'char')
            % dataInName is a file name, can be .mat or a directory

                % dataInName should be a folder with a)matlab files b)tiff files or c) folders 

                    % no matlab files, process individual folders
                    disp('Read images from folders and return matlab data')
                    dir1                        = dir(dataInName);

                    tempValForDotFiles          = 0;
                    handles.numFrames           = size(dir1,1);

                    for counterDir=1:handles.numFrames
                        tempDir=dir1(counterDir).name;
                        if ~strcmp(tempDir(1),'.')
                            dataInName1         =  strcat(dataInName,'/',tempDir);
                            if isdir(dataInName1)
                                %disp(strcat('Processing folder: ',dataInName1));
                                %dataOutName1        =  strcat(dataOutName,'/',tempDir);
                                %restrict to tiff files for the time being
                                dir2                = dir(strcat(dataInName1,'/*.tif'));
                                numSlices           = size(dir2,1);
                                for counterSlice=1:numSlices
                                    tempDir2        = dir2(counterSlice).name;
                                    dataInName2     = strcat(dataInName1,'/',tempDir2);
                                    %if data is an image, read with imread, if not, it may be PIC,
                                    if (strcmp(dataInName2(end-2:end),'PIC'))||(strcmp(dataInName2(end-2:end),'pic'))
                                        % BioRad PIC Format
                                        fid                 = fopen(dataInName2);
                                        [datavect,count]    = fread(fid);
                                        
                                        xDim                = datavect(1)+datavect(2)*256;
                                        yDim                = datavect(3)+datavect(4)*256;
                                        zDim                = datavect(5);

                                        %---- fixed, started reading on element 77
                                        initial=77;
                                        total=76+xDim*yDim*zDim;
                                        dataIn(:,:,counterSlice,counterDir-tempValForDotFiles) = reshape(datavect(initial:total),xDim,yDim,zDim);
                                        
                                    else
                                        dataIn(:,:,counterSlice,counterDir-tempValForDotFiles) = imread(dataInName2);
                                    end
                                end
                                %----- the images read are saved to a file HERE ------
                                %save(dataOutName1,'dataIn');
                                %-----------------------------------------------------
                            else
                                %if data is an image, read with imread, if not, it may be PIC, will read
                                if (strcmp(dataInName1(end-2:end),'PIC'))||(strcmp(dataInName1(end-2:end),'pic'))
                                    % BioRad PIC Format
                                    fid                 = fopen(dataInName1);
                                    [datavect,count]    = fread(fid);
                                    
                                    xDim                = datavect(1)+datavect(2)*256;
                                    yDim                = datavect(3)+datavect(4)*256;
                                    zDim                = datavect(5);
                                    
                                    %---- fixed, started reading on element 77
                                    initial=77;
                                    total=76+xDim*yDim*zDim;
                                    dataIn(:,:,:,counterDir-tempValForDotFiles) = reshape(datavect(initial:total),xDim,yDim,zDim);
                                else
                                    
                                    dataIn(:,:,:,counterDir-tempValForDotFiles) = imread(dataInName1);
                                end
                                
                            end
                        else
                            tempValForDotFiles  = tempValForDotFiles+1;
                        end
                    end
        else
            %----- dataInName is not a char,
            disp('Argument received is not a string, please try again');
            dataIn=[];
            return;
        end
    otherwise
        %----- two arguments received, read files from a folder, and place in a single matlab matrix
        disp('Two arguments received, please try again');
        dataIn=[];
        return;

end
