function [h2d,xhist]=hist2d(data,nbins,toplot)
%-------HIST2D  script that displays the histogram of a 1,2,3D image ----
%-------        reshapes the image: 1d vector and  plots histogram ------
%------------------------------------------------------------------------
%----------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro
%------             PHD     the University of Warwick
%------  Supervisor :   Abhir Bhalerao    -----------
%------  6 November 2001 ----------------------------
%----------------------------------------------------

if nargin<1  help hist2d;  return; end;
%if nargin==1  nbins=50; end;
if ~(isa(data,'double')) data=double(data); end

if nargin==1  nbins=min(min(data(:))):max(max(data(:))); end;
if (nargin==2&nbins==1)  nbins=min(min(data)):max(max(data)); toplot=1; end;
if ~exist('toplot') toplot=0;   end
[lins,cols,levels]=size(data);



%------ revise the case ------
%------    1 1D use hist for either line or column data
%------    2 1D but not line or column, stored in various z
%------    3 2D reshape matrix to vector
%------    4 3D reshape matrix to vector

if ((cols==1|lins==1)&levels==1)
       [h2d,xhist]=hist(data,nbins);
elseif (cols==1&lins==1&levels~=1)
       h2d=permute(data,[2 3 1]);
       xhist=1:levels;
elseif(lins~=1&cols~=1&levels==1)
       data2=reshape(data,1,cols*lins);
       [h2d,xhist]=hist(data2,nbins);
else
       data2=reshape(data,1,cols*lins*levels);
       [h2d,xhist]=hist(data2,nbins);
end;
if (toplot==1)
    bar(xhist,h2d);
    axis tight
end
