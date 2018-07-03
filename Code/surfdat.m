function hand1=surfdat(x,y,datatoplot)
% function surfdat(x,y,datatoplot)
% SURFDAT  Display data received with surf in
% the format desired,
%----------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro
%------             PHD     the University of Warwick
%------  Supervisor :   Abhir Bhalerao    -----------
%------  18 October 2001 ----------------------------
%----------------------------------------------------

%------ no input data is received, error -------------------------
if nargin < 1  help surfdat;  return; end;

%------ Only data to be plotted is received ----------------------
if (nargin==1)|((nargin==2)&(isa('y','char')))|((nargin==3)&(isa('y','char')))
    %h=gca; cla;
    if ~isa(x,'double') x=double(x); end
    [rows,cols,levs]=size(x);
    if (rows~=1&cols~=1&levs~=1)   %----- this is the case of 3D DATA
        if ~isreal(x) x=abs(x); end;
        h=gca;
        %cla;
        hold on;
        plotAxis=0;%x=permute(x,[3 2 1]);
        [rows,cols,levs]=size(x);
        [X,Y]=meshgrid([1:cols],[1:rows]);
        Z=ones(rows,cols);
        %plot3(X,Y,Z,'k.','markersize',1);
        %plot3(Z,Y,X,'k.','markersize',1);
        %plot3(X,Z*cols,Y,'k.','markersize',1);
        thresholdToPlot=max(x(:))/2;
        if ~exist('datatoplot')
             markerS=5;
        else
            markerS=datatoplot;
        end

        %markerS=5;
        if ~exist('y')
            markerCode='ro';markerCode2='g.';
        elseif strcmp(y,'all')
            markerCode = jet;%.^2;
%            markerCode = {0.5*[1 1 1],'bx', 'ro', 'kd', 'c.', 'y+', 'mv', 'g*',...
 %                         'b^', 'r<', 'k>', 'ch', 'yp', 'mo', 'g+'};
                  
            [yhist,xhist]=hist(x(:),[1:max(x(:))]);
            yhistT=find(yhist);
            markerStep=round(linspace(1,64,size(yhistT,2)));
            [changeOrd,indexOrd]=sort(rand(size(markerStep,2),1));
            markerStep=markerStep(indexOrd);
            for cCurrValue=1:numel(yhistT)
                %surfdat(x==cCurrValue,markerCode{1+mod(cCurrValue,14)}); 
                surfdat(x==yhistT(cCurrValue),markerCode(markerStep(cCurrValue),:),markerS); 
            end
            return;
        else
            markerCode=y;markerCode2='g.';
        end
        for cLev=1:levs
            if any(any(x(:,:,cLev)>thresholdToPlot))
                X1=X.*(x(:,:,cLev)>thresholdToPlot); X1(X1==0)=[];
                Y1=Y.*(x(:,:,cLev)>thresholdToPlot); Y1(Y1==0)=[];
                Z1=(levs-cLev+1)*(x(:,:,cLev)>thresholdToPlot); Z1(Z1==0)=[];
                %  plot3(Y1,X1,Z1,'.','color',0.3*[1 1 1 ],'markersize',markerS);
                if isa(markerCode,'char')
                    plot3(Y1,X1,Z1,markerCode,'markersize',markerS);
                else
                    plot3(Y1,X1,Z1,'color',markerCode,'marker','o','markersize',markerS,'linestyle','none');
                end
                %                 X2=X.*(x(:,:,cLev)>thresholdToPlot*1.75); X2(X2==0)=[];
                %                 Y2=Y.*(x(:,:,cLev)>thresholdToPlot*1.75); Y2(Y2==0)=[];
                %                 Z2=(levs-cLev+1)*(x(:,:,cLev)>thresholdToPlot*1.75); Z2(Z2==0)=[];
                %                 plot3(Y2,X2,Z2,markerCode2,'markersize',markerS);

            end
            if exist('Y1')&(plotAxis==1)
                if ~isempty(Y1)
                    plot3(ones(1,levs)*Y1(1),ones(1,levs)*X1(1),[1:levs],'g-.')
                    plot3(ones(1,cols)*Y1(1),[1:cols],ones(1,cols)*Z1(1),'g-.')
                    plot3([1:rows],ones(1,rows)*X1(1),ones(1,rows)*Z1(1),'g-.')
                    plotAxis=0;
                end
            end
        end
        hold off;
        axis([1 rows 1 cols 1 levs]);
        view(85,20);
        set(h,'ZTick',[1:ceil(rows/16):levs])
        set(h,'ZTickLabel',[levs:-1*ceil(rows/16):1])
        rotate3d on;
        grid on;
    elseif xor(rows~=1,cols~=1)&levs==1        %------ This is the case of 1D DATA
        %if ~isreal(x)
        h1=gcf; clf;
        lenX=max(rows,cols);
        x=reshape(x,1,lenX);
        imagX=imag(x);realX=real(x);magX=abs(x);phaX=angle(x);
        maxR=1+max(realX);maxI=1+max(imagX);
        axZeros=zeros(1,lenX);axOnes=ones(1,lenX);axT=1:lenX;
        X=[axZeros;realX];
        Y=[axZeros;imagX];
        Z=[axT;axT];
        h1=gcf;
        subplot(211);     %-----  first plot will include real and imag
        hold off;
        plot3(axT,imagX,realX,'ro','markersize',5);     %----  Main complex function
        hold on;
        plot3(axT,imagX,realX,'r:');                    %---- envolvente
        line(Z,Y,X,'color','b');                        %---- stem of complex
        plot3(axT,axZeros,axZeros,'b--');            %---- axis
        %----- imaginary marginal
        plot3(axT,imagX,-maxR*1.5*axOnes,'m:o');
        %line(Z,[axZeros;imagX],[-maxR*1.5*axOnes;-maxR*1.5*axOnes],'color','c');
        %----- real marginal
        stem3(axT,maxI*1.5*axOnes,realX,'m');
        grid on;axis tight;
        view(15,45);
        %rotate3d on;
        hold off;
        ylabel('imaginary');zlabel('real');
        subplot(212)  %-----  second plot will include magnitude and phase
        hold off;
        plot3(axT,imagX,realX,'ro','markersize',5);     %----  Main complex function
        hold on;
        plot3(axT,imagX,realX,'r:');                    %---- envolvente
        line(Z,Y,X,'color','b');                        %---- stem of complex
        plot3(axT,axZeros,axZeros,'b--');               %---- axis
        %----- phase as marginal
        plot3(axT,phaX*maxI/4,-maxR*10.5*axOnes,'k:o');
        %line(Z,[axZeros;phaX*maxI/4],[-maxR*1.5*axOnes;-maxR*1.5*axOnes],'color','g')
        %----- magnitude as marginal
        stem3(axT,maxI*1.5*axOnes,magX,'k');

        grid on;axis tight
        view(15,45);
        %rotate3d on;
        hold off;
        ylabel('Phase');zlabel('Magnitude')


        %else
        % plot(x);
        %end;
    elseif (rows==1)|(cols==1)|(levs==1)        %------ This is the case of 2D DATA
        x=squeeze(x);

        if ~isreal(x) x=abs(x); end;
        minx=min(min(x));
        maxx=max(max(x));
        data=64*(x-minx)/max(0.01,(maxx-minx));   %----- the maximum ?????
        %         if rows==1&cols~=1&levs~=1
        %             data=permute(data,[2 3 1]);
        %         elseif rows~=1&cols==1&levs~=1
        %             data=permute(data,[1 3 2]);
        %         end
        if nargout>0
            hand1= imagesc(data);
        else
            cla;
            imagesc(data);
        end
        %axis off;
        %axis square;
        zoom on;

    end
    %------- Received X, Y, and data to plot
elseif nargin==3
    if ~isa(datatoplot,'double') datatoplot=double(datatoplot); end
    if ~isreal(datatoplot) datatoplot=abs(datatoplot); end
    image(x,y,datatoplot);
    axis off;
    axis square;
end
%rotate3d on;