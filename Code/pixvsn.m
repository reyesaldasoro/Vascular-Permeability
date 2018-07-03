function [dataout,neighValue,indexValue]=pixvsn(data,nhood,nClass)
% function [dataout,neighValue,indexValue]=pixvsn(data,nhood,nClass)  
% compares value of pixel versus neighbours and change if necessary
%----------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro
%------             PHD     the University of Warwick
%------  Supervisor :   Abhir Bhalerao    -----------
%------  5 November 2001 ----------------------------
%----------------------------------------------------
%------ In case that more than 6 neighbours have a different
%------ value, set to the value of the neighbours
%------ This is done by shifting to all directions and comparing with data
%------ Input  : Data with neighbours to compare
%------          sqSize is used in case data is in squares (optional) 
%------ Output : Data with uniformed neighbours
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------

if nargin<1;                help pixvsn;        return;     end
if ~(isa(data,'double'));   data=double(data);              end
minData                     = min(data(:));
if minData<=0;              data = data-minData+1;          end
[rows,cols,levels]          = size(data);
methodSelection             = 3;
%-----------------------------------------------------------------
%------ define sizes and neighbourhoods if needed, normally 1 ----
%-----------------------------------------------------------------
if nargin==1    
    if levels==1;           nhood=[3 3];    % nClass=4;  %------ 2D data
    else                    nhood=[3 3 3];  % nClass=4;  %------ 3D data    %else nhood=[sqSize+2 sqSize+2]; 
    end 
    [y,x]                   = hist2d(data);
    nClass                  = x(end);
end
if nargin==2
    if all(size(nhood)==1)          %------ nClass given
        nClass              = nhood;
        if levels==1;       nhood=[3 3]; else nhood=[3 3 3]; end
    else                            %------ nhood given
        [y,x]               = hist2d(data);
        nClass              = x(end);        
    end
end
%---------------------------------------------------------------------------------
%------ it can be performed by rearrangement into vol to col or block process ----
%---------------------------------------------------------------------------------
if xor(levels==1,length(nhood)==2)
    disp('Data and neighbourhood are not consistent in dimensions');
    dataout=[];neighValue=[];%indexValuereturn=[];
    return;
end

%------ Pad the region with zeros -----------------
padData = zeros(size(data) +  (nhood-1));
if levels ==1
    padData(floor((nhood(1)-1)/2)+(1:rows),floor((nhood(2)-1)/2)+(1:cols)) = data;
else
    padData(floor((nhood(1)-1)/2)+(1:rows),floor((nhood(2)-1)/2)+(1:cols),floor((nhood(3)-1)/2)+(1:levels)) = data;
end

if methodSelection==1
    %---------------------------------------------------------------------------------
    %-----------------    2D or 3D re-arranging the matrix   -------------------------
    %---------------------------------------------------------------------------------
    block=nhood;
    ma=rows; na=cols;pa=levels;         %[ma,na,pa] = size(data);
    if (length(block)<3); block(3)=1; end
    
%     %------ Pad the region with zeros -----------------
%     padData = zeros(size(data) +  (nhood-1));
%     padData(floor((block(1)-1)/2)+(1:ma),floor((block(2)-1)/2)+(1:na),floor((block(3)-1)/2)+(1:pa)) = data;
%     %------ Convert neighborhoods of matrix A to columns -----
    x = vol2col(padData,nhood);
    %------ Form weight vector, from there to determine isolated/blob/edge/region
    if pa==1                                        %------ 2D data
        %weightVector=[1 2 1 2 2.5 2 1 2 1]';
        weightVector=[1 1 1 1 1 1 1 1 1]';
    else                                            %------ 3D data
        weightVector=[ones(4,1);2;ones(5,1);2;1;2;2.5;2;1;2;ones(5,1);2;ones(4,1)];
    end
    weightMat=repmat(weightVector,[1 size(x,2)]);
    %------ Classification process ---------------------------  
    for cClass=1:nClass
        numNeighInClass(cClass,:)=sum((x==cClass).*weightMat);
        %numNeighInClass(cClass,:)=sum((x==cClass));
    end
    [maxclass,class]=max(numNeighInClass);
    dataout=reshape(class,[ma na pa]);
    %------ in case they are requested the index values can be provided
    if nargout>1
        indexValue=zeros(ma,na,pa,nClass);
        if pa==1
            for cClass=1:nClass
                indexValue(:,:,cClass)=reshape((numNeighInClass(cClass,:)),[ma na pa]);
            end
        else
            for cClass=1:nClass
                indexValue(:,:,:,cClass)=reshape((numNeighInClass(cClass,:)),[ma na pa]);
            end
        end
   end
elseif methodSelection==2
    clear numNeighInClass
    %---------------------------------------------------------------------------------
    %-----------------    block process                      -------------------------
    %---------------------------------------------------------------------------------
    
    weightVector                        = [ones(4,1);2;ones(5,1);2;1;2;2.5;2;1;2;ones(5,1);2;ones(4,1)];
    weightMat                           = reshape(weightVector,3,3,3);
    dataout                             = zeros(size(data));
    numNeighInClass(nClass)             = 0;
    for cRows=1:rows
        for cCols=1:cols
            for cLevs=1:levels
                blockData=reshape(padData(cRows:cRows+nhood(1)-1,cCols:cCols+nhood(2)-1,cLevs:cLevs+nhood(3)-1),1,27);
                for cClass=1:nClass
                    %numNeighInClass(cClass)=sum(sum(sum((blockData==cClass).*weightMat)));
                    numNeighInClass(cClass)=((blockData==cClass)*weightVector);
                end
                [maxCl,dataout(cRows,cCols,cLevs)]=max(numNeighInClass);
            end
        end
    end
else
    %---------------------------------------------------------------------------------
    %-----------------    whole block process                      -------------------------
    %---------------------------------------------------------------------------------
    %weightVector=[ones(4,1);2;ones(5,1);2;1;2;2.5;2;1;2;ones(5,1);2;ones(4,1)];  %----weigth vector for weighted comparison
    weightVector=[ones(4,1);4;ones(5,1);4;1;4;8.5;4;1;4;ones(5,1);4;ones(4,1)];  %----stronger weight vector for weighted comparison
    %weightVector=[ones(27,1)];                                                  %----weigth vector for non weighted comparison
    %weightVector=[ones(9,1);zeros(9,1);ones(9,1)];                              %----weigth vector for interslice comparison
    weightMat=reshape(weightVector,3,3,3);
        indexValue=zeros(rows,cols,levels,nClass);

    if levels==1
        for cRows=-1:1
            for cCols=-1:1
                blockData=padData(cRows+2:cRows+rows+1,cCols+2:cCols+cols+1);
                for cClass=1:nClass
                    %indexValue(:,:,:,cClass)=indexValue(:,:,:,cClass)+(blockData==cClass)*weightMat(cRows+2,cCols+2,2);
                    indexValue(:,:,:,cClass)=indexValue(:,:,:,cClass)+(blockData==cClass);
                end
                
            end
        end
    else                    %------ enter here only for 3D Cases
        
        for cRows=-1:1
            for cCols=-1:1
                for cLevs=-1:1
                    blockData=padData(cRows+2:cRows+rows+1,cCols+2:cCols+cols+1,cLevs+2:cLevs+levels+1);
                    for cClass=1:nClass
                        indexValue(:,:,:,cClass)=indexValue(:,:,:,cClass)+(blockData==cClass)*weightMat(cRows+2,cCols+2,cLevs+2);
                        %indexValue(:,:,:,cClass)=indexValue(:,:,:,cClass)+(blockData==cClass);
                    end
                    
                end
            end
        end
    end
        [maxCl,dataout]=max(indexValue,[],4);
    
%sum(sum(sum(dataout~=dataout1)))
    
end

%---- neighValue keeps the number of neighbour pixels that are same as itself
%---- now witha loop but think of something more clever!!!!!!
%---- 11 August 2003 developed for the Pyramid classification of Spann/Wilson
neighValue(rows,cols,levels)=0;
if nargout>=2
    for cRows=1:rows
        for cCols=1:cols
            for cLevs=1:levels
                neighValue(cRows,cCols,cLevs)=indexValue(cRows,cCols,cLevs,data(cRows,cCols,cLevs));
            end
        end
    end
end
        


if minData<=0; dataout=dataout+minData-1; end
