% --- calculates the classical sleep metrics for a population of flies over
%     the duration of an experiment (sleep bouts, sleep duration and
%     average sleep bout duration per hour)
function pData = SleepMetric(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Population Sleep Metrics';
pData.Type = {'Pop','Multi'};
pData.fType = [1 1 1 1];
pData.rI = initFuncReqInfo(pData);

% initialises the other fields  (if input argument provided)
if (nargin == 1)    
    % parameter data struct initialisation
    snTotL = snTot(1);
    pData.cP = initCalcPara(snTot);
    pData.pP = initPlotPara(snTot);
    pData.oP = initOutputPara(snTot);
    pData.pF = initPlotFormat(snTotL);
    
    % sets the apparatus name/count
    pData.appName = snTotL.iMov.pInfo.gName;
    pData.nApp = length(pData.appName);   
    
    % special parameter fields/data struct
    isLong = ~detIfShortExpt(field2cell(snTotL,'T'));
    [pData.hasSP,pData.hasRS,pData.hasRC] = deal(true,false,isLong);
    pData.sP = initSpecialPara(snTotL,pData); 
end
    
% ----------------------------------------------------------------------- %
% ---                 PARAMETER STRUCT SETUP FUNCTIONS                --- %
% ----------------------------------------------------------------------- %

% --- sets the function required information struct
function rI = initFuncReqInfo(pData)

% memory allocation
rI = struct('Scope',[],'Dur',[],'Shape',[],...
            'Stim',[],'Spec',[],'SpecFcn',[],'ClassicFcn',true);

% sets the struct fields
rI.Scope = setFuncScopeString(pData.Type);
rI.Dur = 'None';
rI.Shape = 'None';
rI.Stim = 'None';
rI.Spec = 'None';        
        
% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 1;
sz0 = [1 2 4 6 8 12 24];
cP = setParaFields(nPara);

% determines if the experiment is long (if not reduce the group size list)
isLong = ~detIfShortExpt(field2cell(snTot,'T'));
if ~isLong; sz0 = sz0(1); end

% sets the tab list names
a = {'1 - General'};

% sets the parameter list for 
pList = cellfun(@num2str,num2cell(sz0),'un',0);

% sets the parameter fields
cP(1) = setParaFields(a{1},'List',{1+isLong,pList},'nGrp',...
    'Number of Daily Time Groups',[],{1,~isLong});
cP(1).TTstr = 'The number of groups that the day is split up into';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 7;                       
pP = setParaFields(nPara);
isLong = ~detIfShortExpt(field2cell(snTot,'T'));

% sets the parameter list strings
pList = {'Sleep Bouts','Sleep Duration','Mean Bout Duration','Combined Bout & Duration'};
pList2 = {'Bar Graph','Boxplot'};
pList3 = {'Hourly','Half-Daily','Daily'};

% sets the tab list names
a = {'1 - General','2 - Error/Outliers'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},'Number',0.75,'pW','Bar Plot Relative Width',[0 1 false],{3,1});
pP(2) = setParaFields(a{1},'List',{1,pList},'pMet','Plot Metrics');
pP(3) = setParaFields(a{1},'List',{1,pList2},'pType','Plot Type');
pP(4) = setParaFields(a{1},'List',{1,pList3},'kRate','Sleep Bout Equivalent Rate',[],{2,[1 4]});
pP(5) = setParaFields(a{1},'Boolean',isLong,'pltDN','Show Day/Night Background',[],{0,~isLong});
pP(6) = setParaFields(a{1},'Boolean',0,'plotGrid','Show Axis Gridlines');
pP(7) = setParaFields(a{2},'Boolean',1,'plotErr','Show Error Bars/Outliers');

% sets the tool-tip strings
pP(4).TTstr = 'The rate over which the sleep bouts are calculated';

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
nApp = length(snTot.iMov.ok);    
pF = setFormatFields(nApp);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],'',1);
pF.yLabel = setFormatFields([],'',1);
pF.zLabel = setFormatFields([],'Sleep (Minutes/Hour)',1);
pF.Axis = setFormatFields([],[]);

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Title(i).String = snTot.iMov.pInfo.gName{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);
[Stats,Type] = deal({'CompMulti','Tgrp'},[1 3]);

% sets the independent variable fields
oP = addXVarField(oP,'Time Group','Tgrp','Group');

% sets the dependent variable fields
oP = addYVarField(oP,'Sleep Bout','N',Stats,Type,{'Tgrp'},1);
oP = addYVarField(oP,'Sleep Duration','Ts',Stats,Type,{'Tgrp'},1);
oP = addYVarField(oP,'Bout Duration','R',Stats,Type,{'Tgrp'},1);

% ----------------------------------------------------------------------- %
% ---                       CALCULATION FUNCTION                      --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [plotD,ok] = calcFunc(snTot,pData,gPara,cP,varargin)

% initialises the calculation parameters (if not already initialised)
if (nargin == 3)
    % retrieves the parameter struct
    cP = retParaStruct(pData.cP,gPara);    
end

% sets the movement type (based on the global parameters)
if (strcmp(gPara.movType,'Absolute Location'))
    cP.movType = 'Absolute Range';
end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);
nGrp = str2double(cP.nGrp);

% fixed parameters
cP.tBin = 60;
a = NaN(nGrp,1);

% sets the daily time group strings
Tgrp = setTimeGroupStrings(nGrp,cP.Tgrp0);                     

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,NaN(1,nGrp),...
                                 'N',[],'N_mn',a,'N_sem',a,...
                                 'Ts',[],'Ts_mn',a,'Ts_sem',a,...
                                 'R',[],'R_mn',a,'R_sem',a,...                                    
                                 'Tgrp',Tgrp);                                                                                
                             
% creates the waitbar figure
wStr = {'Overall Progress','Sleep Metric Calculations'};
wOfs = (nExp > 1); wStr = wStr((2-wOfs):end);
if nargin == 5
    % resets the strings
    [h,wOfs] = deal(varargin{1},wOfs+1);
    
    % resets the progressbar fields
    wStr = [h.wStr(1),wStr];
    for i = 2:length(wStr); h.Update(i,wStr{i},0); end
else
    % creates the waitbar figure
    h = ProgBar(wStr,'Sleep Metric Calculations');    
end

% --------------------------- %
% --- METRIC CALCULATIONS --- %
% --------------------------- %

% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure (if more than one solution file)
    if nExp > 1
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr{wOfs},i,nExp);
        if h.Update(wOfs,wStrNw,i/(1+nExp))
            [plotD,ok] = deal([],false);
            return
        else
            h.Update(1+wOfs,wStr{1+wOfs},0);
        end
    end

    % sets the total time and ok flags
    Ttot = cell2mat(snTot(i).T);
    flyok = snTot(i).iMov.flyok;
    
    % determines the binned indices (for time length, tBin) and determines 
    % the bins which has at least two time points
    h.Update(1+wOfs,'Determining Time Bins',0.50);
    indB = detTimeBinIndices(Ttot,cP.tBin);
    indB(cellfun(@length,indB) < cP.tBin/2) = {[]};        
    
    % calculates the sleep metrics
    for j = 1:nApp
        % updates the waitbar figure
        wStrNw = sprintf(['Calculating Sleep Metrics (Region ',...
                          '%i of %i)'],j,nApp);
        if h.Update(1+wOfs,wStrNw,0.50*(1+(j/(nApp+1))))
            % if the user cancelled, then exit the function
            [plotD,ok] = deal([],false);
            return            
        end
                 
        % only calculate if values exist...
        if ~isempty(snTot(i).Px{j})      
            % calculates the sleep metrics
            [nBoutNw,tSleepNw] = ...
                    calcSleepMetrics(snTot(i),Ttot,indB,cP,j,flyok{j});
            RNw = cellfun(@(x,y)(x./y),tSleepNw,nBoutNw,'un',0);
            
            % sets the raw values into the storage arrays
            plotD(j) = setRawDataValues(...
                            plotD(j),snTot(i),nBoutNw,[],'N',i,j,2);
            plotD(j) = setRawDataValues(...
                            plotD(j),snTot(i),tSleepNw,[],'Ts',i,j,2);
            plotD(j) = setRawDataValues(...
                            plotD(j),snTot(i),RNw,[],'R',i,j,2);                
        end
    end
    
    % updates the waitbar figure
    h.Update(1+wOfs,'Sleep Metric Calculation Complete!',1);    
end
   
% calculates the metrics for all apparatus
for j = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf('Calculating Sleep Metric (Region %i of %i)',j,nApp);
    if h.Update(1+wOfs,wStrNw,(j+1)/(2+nApp))
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    end               
        
    % calculates the metric statistics for each of the values
    plotD(j) = calcMetricStats(plotD(j),'N');
    plotD(j) = calcMetricStats(plotD(j),'Ts');
    plotD(j) = calcMetricStats(plotD(j),'R'); 
end

% closes the waitbar figure
if ~h.Update(1+wOfs,'Sleep Metric Calculations Complete!',1)
    if nargin < 5; h.closeProgBar(); end
end

% ----------------------------------------------------------------------- %
% ---                        PLOTTING FUNCTION                        --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,varargout] = plotFunc(snTot,pData,plotD,iPlot)

% retrieves the plotting paraeter struct
pP = retParaStruct(pData.pP);
sP = retParaStruct(pData.sP);
cP = retParaStruct(pData.cP);
pF = pData.pF;

% if the incorrect combination is used, then exit with an error (not
% possible to plot the day/night separation and combined histograms
% together)
if (strcmp(pP.pMet,'Combined Bout & Duration') && strcmp(pP.pType,'Boxplot'))
    eStr = 'Not possible to use boxplot with Combined Bout & Duration.';
    waitfor(errordlg(eStr,'Incorrect Plot Type','modal'))
    return
end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the plotting indices and subplot indices
[ind,m,n] = deal(find(sP.Sub.isPlot),sP.Sub.nRow,sP.Sub.nCol);
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% sets the time group value and the x-plot values
nGrp = str2double(cP.nGrp);
xi = 1:nGrp;

% other parameters
[sMlt,nTick,fAlpha] = deal(1,5,0.4);
[ix,iy,tRot] = deal([1 1 2 2],[1 2 2 1],30);
isComb = strcmp(pP.pMet,'Combined Bout & Duration');

% sets the output arguments
varargout{1} = pData;

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% retrieves the formatting struct
if (nGrp == 1)
    pF = retFormatStruct(pF,1);
else
    if (isempty(m)); szMx = 1; else szMx = max([m n]); end
    pF = retFormatStruct(pF,szMx);
end

% sets the y-label string
if (nargin == 4)
    switch (iPlot)
        case (1)
            switch (pP.kRate)
                case ('Hourly')
                    pF.yLabel.String = sprintf('Sleep (Bouts/Hour)');
                case ('Half-Daily')
                    sMlt = 12;
                    pF.yLabel.String = sprintf('Sleep (Bouts/%sDay)',char(189));
                case ('Daily')
                    sMlt = 24;
                    pF.yLabel.String = sprintf('Sleep (Bouts/Day)');
            end
        case (2)
            pF.yLabel.String = sprintf('Sleep (Minutes/Hour)');            
        case (3)
            pF.yLabel.String = sprintf('Sleep (Minutes/Bout)');
    end
else
    switch (pP.pMet)
        case {'Sleep Bouts','Combined Bout & Duration'}
            switch (pP.kRate)
                case ('Hourly')
                    pF.yLabel.String = sprintf('Sleep (Bouts/Hour)');
                case ('Half-Daily')
                    sMlt = 12;
                    pF.yLabel.String = sprintf('Sleep (Bouts/%sDay)',char(189));
                case ('Daily')
                    sMlt = 24;
                    pF.yLabel.String = sprintf('Sleep (Bouts/Day)');
            end            
        case ('Sleep Duration')
            pF.yLabel.String = sprintf('Sleep (Minutes/Hour)');
        case ('Mean Bout Duration')
            pF.yLabel.String = sprintf('Sleep (Minutes/Bout)');
    end
end

% removes the label indices
[pF.xLabel.ind,pF.yLabel.ind,pF.zLabel.ind] = deal(NaN);
    
% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% initialisations
hAxB = cell(nApp,1);
xLim = [1 nGrp]+0.5*[-1 1];

% sets the fill coordinates for the day/night comparison
yLim = 0;
[xFillD,xFillN,yFill] = deal([xLim(1) mean(xLim)],[mean(xLim) xLim(2)],[0 1]);
[hErr,hBar] = deal(cell(nApp,1));

% sets the plot 
switch (pP.pMet)
    case ('Sleep Bouts') % case is the sleep bout count
        [yL,pStr,nTick] = deal([],'N',6);        
    case ('Sleep Duration') % case is the sleep duration
        [yL,pStr] = deal([0,60],'Ts');
    case ('Mean Bout Duration') % case is the avg bout duration
        [yL,pStr] = deal([],'R');        
end    

% loops through all the indices plotting the solutions
if (nGrp == 1)
    [xi,pF.Title(1).String] = deal(1:nApp,pP.pMet);
    [xLim,m,n] = deal(xi([1 end])+0.5*[-1 1],1,1); 
    
    % sets up the plot axis
    hAxB = createSubPlotAxes(hP);    
    axis(hAxB,'on'); hold on; 
    
    % sets the mean/SEM metric values
    if strcmp(pP.pMet,'Combined Bout & Duration')
        % plots the double axis bar
        [N,T] = deal(field2cell(p,'N_mn',1),field2cell(p,'Ts_mn',1));
        [Ns,Ts] = deal(field2cell(p,'N_sem',1),field2cell(p,'Ts_sem',1));
        
        if (pP.plotErr)
            [hAxB,hBar] = plotDoubleAxisBar(hAxB,xi,sMlt*N,...
                        T,pP.pW,[],nTick,[NaN 60],sMlt*Ns,Ts);                             
        else
            [hAxB,hBar] = plotDoubleAxisBar(hAxB,xi,sMlt*N,...
                        T,pP.pW,[],nTick,[NaN 60],[],[]);                             
        end
    else
        % creates the bar graph/boxplot
        plotBarBoxMetrics(hAxB,xi,p,pStr,pP,yL,'b',sMlt);    
    end
    
    % updates the axis properties
    set(hAxB,'xticklabel',[],'xlim',xLim,'linewidth',1.5,'UserData',1);
    if (isComb)
        % formats the double plot axis
        pF = formatDoubleAxis(hAxB,hBar,pF,1);     
    else
        % formats the single axis
        formatPlotAxis(hAxB,pF,1);                                 
    end              

    % turns the grid on (if specified)
    if (pP.plotGrid); set(hAxB(1),'ygrid','on'); end     
else
    for j = 1:nApp
        % initialises the subplot axes
        [i,hAxB{j}] = deal(ind(j),createSubPlotAxes(hP,[m,n],j));
        hold(hAxB{j},'on');
        
        % plots the background image
        if (pP.pltDN)
            patch(xFillD(ix),yFill(iy),'y','FaceAlpha',fAlpha,'tag','hDN')
            patch(xFillN(ix),yFill(iy),'k','FaceAlpha',fAlpha,'tag','hDN')        
        end

        % sets the mean/SEM metric values
        if (isComb)
            % plots the double axis bar
            if (pP.plotErr)
                [hAxB{j},hBar{j},hErr{j}] = plotDoubleAxisBar(hAxB{j},xi,sMlt*p(j).N_mn,...
                            p(j).Ts_mn,pP.pW,[],nTick,[NaN 60],...
                            sMlt*p(j).N_sem,p(j).Ts_sem);                             
            else
                [hAxB{j},hBar{j},hErr{j}] = plotDoubleAxisBar(hAxB{j},xi,sMlt*p(j).N_mn,...
                            p(j).Ts_mn,pP.pW,[],nTick,[NaN 60],[],[]);                             
            end
        else
            % creates the bar graph/boxplot
            yLimNw = plotBarBoxMetrics(hAxB{j},xi,p(j),pStr,pP,yL,'b',sMlt);
            if (any(strcmp(pP.pMet,{'Mean Bout Duration','Sleep Bouts'})))
                yLim = max(yLim,yLimNw);                
            else
                yLim = [0 60];
            end
        end

        % updates the axis properties
        set(hAxB{j},'xticklabel',[],'xlim',xLim,'linewidth',0.5);
        if (isComb)            
            % formats the double plot axis
            pF = formatDoubleAxis(hAxB{j},hBar{j},pF,i);                 
        else
            % formats the single axis
            formatPlotAxis(hAxB{j}(1),pF,i);                             
        end      

        % turns the grid on (if specified)
        if (pP.plotGrid); set(hAxB{j}(1),'ygrid','on'); end           
    end
end
        
% ensures that the D/N 
if (nGrp > 1)
    if (~isempty(hAxB))
        if (strcmp(pP.pMet,'Combined Bout & Duration'))
            % determines the axis with the largest values 
            yLimMx = cellfun(@(x)(max(get(x(1),'ylim'))),hAxB);
            [YmxNw,imx] = max(yLimMx);

            % sets the new limits and tick marks
            yLimNw = get(hAxB{imx}(1),'ylim');
            yTickNw = get(hAxB{imx}(1),'ytick');

            % loops through all of the other graphs 
            for i = 1:length(hAxB)
                if (yLimMx(i) ~= YmxNw)
                    yScl = YmxNw/yLimMx(i);
                    set(hAxB{i}(1),'ylim',yLimNw,'ytick',yTickNw)

                    % updates the bar graph y-values
                    hBar2 = findall(hAxB{i}(1),'facecolor','red','type','hggroup');    
                    set(hBar2,'ydata',get(hBar2,'ydata')*yScl);

                    % updates the error bars y-values
                    if (pP.plotErr && ~isnan(hErr{i}(2)))
                        hC = get(hErr{i}(2),'Children');
                        for j = 1:length(hC)
                            set(hC(j),'ydata',get(hC(j),'ydata')*yScl)
                        end
                    end
                end
            end
        else
            if (isempty(yL))
                YmxNw = [yLim(1),detOverallLimit(yLim)];
                if (isnan(YmxNw)); YmxNw = [0,1]; end     
            else
                YmxNw = yL;
            end

            for i = 1:length(hAxB)
                set(hAxB{i},'ylim',YmxNw)
            end
        end

        if (pP.pltDN)
            % sets the upper y-axis limit
            if strcmp(pP.pMet,'Sleep Duration')
                % case is the sleep duration
                YmxNw = 60;                
            else
                % case is ther other metrics
                YmxNw = max(cellfun(@(x)(setStandardYAxis(x,[],nTick,YmxNw)),hAxB));
            end
            
            % retrieves the fill handles
            hFill = cellfun(@(x)(findobj(x(1),'tag','hDN')),hAxB,'un',0);
            yDataNw = [0 YmxNw*[1 1] 0];
            cellfun(@(x)(set(x,'yData',yDataNw)),hFill);

            if (length(hAxB{1}) == 2)
                hFill = cellfun(@(x)(findobj(x(2),'tag','hDN')),hAxB,'un',0);
                cellfun(@delete,hFill)
            end
        end    
    else    
        if (strcmp(pP.pMet,'Sleep Bouts'))
            hAxB = findall(get(hAxB,'parent'),'type','axes');
            for i = 1:length(hAxB)
                yLimNw = setStandardYAxis(hAxB(i),[],nTick,yLim);
            end
        else
            yLimNw = 60;
            set(findall(get(hAxB,'parent'),'type','axes'),'yLim',[0 yLimNw])
        end

        % retrieves the fill handles
        hFill = findall(get(hAxB(1),'parent'),'tag','hDN');
        set(hFill,'yData',[0 yLimNw*[1 1] 0]);     
    end

    % prevents an error with the figure creation
    if (strcmp(pP.pMet,'Sleep Bouts'))    
        figTmp = findall(0,'type','figure','tag','nwPlot');
        if ~isempty(figTmp)
            setObjVisibility(figTmp,'on'); 
            pause(0.05); 
            setObjVisibility(figTmp,'off'); 
        end
    end
end

% sets the non-aligned x/y labels
formatMultiXYLabels(hAxB,pF,[m,n]);  

% resets the axis positions
if (nGrp == 1)
    % case is for a single group    
    resetAxesPos({hAxB},m,n);
        
    % sets the group strings
    grpName = snTot(1).iMov.pInfo.gName(ind);
    setGroupString(hAxB(1),pF(1),xi,grpName,tRot,-0.005);
    
    % resets the 2nd-axis location (if it exists)
    if (length(hAxB) > 1)
        set(hAxB(2),'position',get(hAxB(1),'position'));
    end    
else
    % case is for multiple subplots
    resetAxesPos(hAxB,m,n);
end
    
% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% retrieves the calculation parameters
cP = retParaStruct(pData.cP);

% if only groups, then re-assign the grouping string
if (str2double(cP.nGrp) == 1)
    pData.oP.xVar(1).Type = 'Other';
end