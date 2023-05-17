% --- calculates the waking sleep metrics (waking distance travelled and 
%     waking duration per hour)
function pData = WakingMetric(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Population Waking Metrics';
pData.Type = {'Pop','Multi'};
pData.fType = [1 1 1 1];
pData.rI = initFuncReqInfo(pData);
pData.dcFunc = @dataCursorFunc;

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
    [pData.hasSP,pData.hasRC] = deal(true,isLong);
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
    'Number of Daily Time Groups',[],{0,~isLong});
cP(1).TTstr = 'The number of groups that the day is split up into';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 6;                       
pP = setParaFields(nPara);
isLong = ~detIfShortExpt(field2cell(snTot,'T'));

% sets the parameter list for 
pList = {'Average Total Speed','Average Wake Speed','Waking Duration'};
pList2 = {'Bar Graph','Boxplot'};

% sets the tab list names
a = {'1 - General','2 - Error/Outliers'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},'List',{1,pList},'pMet','Plot Metrics');
pP(2) = setParaFields(a{1},'List',{1,pList2},'pType','Plot Type');
pP(3) = setParaFields(a{1},'Number',0.75,'pW','Bar Plot Relative Width',[0 1 false],{1,1});
pP(4) = setParaFields(a{1},'Boolean',isLong,'pltDN','Show Day/Night Background',[],{0,~isLong});
pP(5) = setParaFields(a{1},'Boolean',0,'plotGrid','Show Axis Gridlines');
pP(6) = setParaFields(a{2},'Boolean',1,'plotErr','Show Error Bars/Outliers');

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
nApp = length(snTot.iMov.ok);    
pF = setFormatFields(nApp);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],'',1);
pF.yLabel = setFormatFields([],'',1);
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
oP = addYVarField(oP,'Total Speed','Dt',Stats,Type,{'Tgrp'},1);
oP = addYVarField(oP,'Wake Speed','Dw',Stats,Type,{'Tgrp'},1);
oP = addYVarField(oP,'Wake Duration','Tw',Stats,Type,{'Tgrp'},1);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(~,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
cP = retParaStruct(dcObj.pData.cP);
pP = retParaStruct(dcObj.pData.pP);
sP = retParaStruct(dcObj.pData.sP);

% sets the common class fields
dcObj.pType = pP.pType;
dcObj.yName = pP.pMet;
dcObj.xName = 'Time Group';
dcObj.xGrp = dcObj.plotD{1}(1).Tgrp;
dcObj.combFig =  strcmp(cP.nGrp,'1');
[dcObj.xUnits,dcObj.yGrp] = deal([]);
dcObj.grpName = dcObj.pData.appName(sP.Sub.isPlot);

% sets the metric units (based on the plot metric)
switch dcObj.pType
    case {'Average Total Speed', 'Average Wake Speed'}
        % case is the average total speed
        dcObj.yUnits = 'mm/sec';
        
    case 'Waking Duration'
        % case is the waking duation
        dcObj.yUnits = 'sec';
        
end

% sets up the data cursor string
dTxt = dcObj.setupCursorString();

% ----------------------------------------------------------------------- %
% ---                       CALCULATION FUNCTION                      --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [plotD,ok] = calcFunc(snTot,pData,gPara,cP,varargin)

% initialises the calculation parameters (if not already initialised)
if nargin == 3
    cP = retParaStruct(pData.cP,gPara);  
end

% sets the movement type (based on the global parameters)
if strcmp(gPara.movType,'Absolute Location')
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
                             'Dt',[],'Dt_mn',a,'Dt_sem',a,...
                             'Dw',[],'Dw_mn',a,'Dw_sem',a,...
                             'Tw',[],'Tw_mn',a,'Tw_sem',a,...
                             'Tgrp',Tgrp,'movType',cP.movType);                                           
                             
% creates the waitbar figure
wStr = {'Overall Progress','Waking Metric Calculations'};
wOfs = (nExp > 1); wStr = wStr((2-wOfs):end);

% sets up/resets the progressbar
if nargin == 5
    % retrieves the input objects
    [h,wOfs] = deal(varargin{1},wOfs+1);
    wStr = [h.wStr(1),wStr];
    
    % resets the waitbar figure
    for i = 2:length(wStr); h.Update(i,wStr{i},0); end
else
    % creates the waitbar figure
    h = ProgBar(wStr,'Wake Metric Calculations');
end

% sets the time multiplier
if strcmp(cP.movType,'Absolute Range')
    mlt = 1/60;
else
    mlt = 1;
end

% --------------------------- %
% --- METRIC CALCULATIONS --- %
% --------------------------- %

% memory allocation
nDayMx = 0;

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
    nDayMx = max(nDayMx,ceil(convertTime(Ttot(end),'sec','days')));

    % determines the binned indices (for time length, tBin) and determines the
    % bins which has at least two time points
    h.Update(1+wOfs,'Determining Time Bins',0.50);
    indB = detTimeBinIndices(Ttot,cP.tBin);
    indB(cellfun('length',indB) < cP.tBin/2) = {[]};
    
    % calculates the sleep metrics
    for j = 1:nApp
        % updates the waitbar figure
        wStrNw = sprintf(['Calculating Waking Metrics (Region ',...
                          '%i of %i)'],j,nApp);
        if h.Update(1+wOfs,wStrNw,0.50*(1+(j/(nApp+1))))
            % if the user cancelled, then exit the function
            [plotD,ok] = deal([],false);
            return                        
        end
         
        % only calculate if values exist...
        if (~isempty(snTot(i).Px{j}))        
            % calculates the wake metrics
            [dTotNw,dWakeNw,tWakeNw] = ...
                    calcWakingMetrics(snTot(i),Ttot,indB,cP,j,flyok{j}); 

            % sets the raw values into the storage arrays
            plotD(j) = setRawDataValues(...
                        plotD(j),snTot(i),dTotNw,[],'Dt',i,j,2,mlt);
            plotD(j) = setRawDataValues(...
                        plotD(j),snTot(i),dWakeNw,[],'Dw',i,j,2,mlt);
            plotD(j) = setRawDataValues(...
                        plotD(j),snTot(i),tWakeNw,[],'Tw',i,j,2);                
        end
    end
    
    % updates the waitbar figure
    h.Update(1+wOfs,'Waking Metric Calculation Complete!',1);    
end
    
% calculates the metrics for all apparatus
for j = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf(['Calculating Waking Metric (Region ',...
                      '%i of %i)'],j,nApp);
    if h.Update(1+wOfs,wStrNw,(j+1)/(2+nApp))
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    end               
             
    % calculates the metric statistics for each of the values
    plotD(j) = calcMetricStats(plotD(j),'Dt');
    plotD(j) = calcMetricStats(plotD(j),'Dw');
    plotD(j) = calcMetricStats(plotD(j),'Tw'); 
end

% closes the waitbar figure
if ~h.Update(1+wOfs,'Waking Metric Calculations Complete!',1)
    if (nargin < 5); h.closeProgBar(); end
end

% ----------------------------------------------------------------------- %
% ---                        PLOTTING FUNCTION                        --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,varargout] = plotFunc(snTot,pData,plotD,ind)

% retrieves the plotting paraeter struct
pP = retParaStruct(pData.pP);
sP = retParaStruct(pData.sP);
cP = retParaStruct(pData.cP);
pF = pData.pF;

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the plotting indices and subplot indices
[ind,m,n] = deal(find(sP.Sub.isPlot),sP.Sub.nRow,sP.Sub.nCol);
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% calculates the upper y-axis limits
nGrp = str2double(cP.nGrp);

% other parameters
[nTick,fAlpha,xi] = deal(5,0.4,1:nGrp);
[ix,iy] = deal([1 1 2 2],[1 2 2 1]);

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% retrieves the formatting struct
if nGrp == 1
    pF = retFormatStruct(pF,1);
else
    if isempty(m); szMx = 1; else; szMx = max([m n]); end
    pF = retFormatStruct(pF,szMx);
end

% sets the y-axis label string
switch pP.pMet
    case 'Waking Duration'
        pF.yLabel.String = 'Wake Duration (min/hour)';
        
    case 'Average Total Speed'
        pF.yLabel.String = 'Average Speed (mm/sec)';
       
    case 'Average Wake Speed'
        pF.yLabel.String = 'Wake Activity (mm/sec)';

end

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');
    
% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% initialisations
hAxB = cell(1,nApp);
xLim = [1 nGrp]+0.5*[-1 1];
[yLim,tRot] = deal(0,30);

% sets the values to be plotted
switch (pP.pMet)
    case ('Average Total Speed') % case is the average total speed
        [yL,pStr] = deal([],'Dt');
    case ('Average Wake Speed') % case is the average wake speed
        [yL,pStr] = deal([],'Dw');
    case ('Waking Duration') % case is the wake duration
        [yL,pStr] = deal([0 60],'Tw');
end 

% day/night background fill object coordinates
[xFillD,xFillN,yFill] = deal([xLim(1) mean(xLim)],[mean(xLim) xLim(2)],[0 1]);

% loops through all the indices plotting the solutions
if (nGrp == 1)
    % re-initialistions
    [xi,pF.Title(1).String] = deal(1:nApp,pP.pMet);
    [xLim,m,n] = deal(xi([1 end])+0.5*[-1 1],1,1); 
    
    % sets up the plot axis
    hAxB = createSubPlotAxes(hP);
    axis(hAxB,'on'); hold on; 
    
    % creates the bar graph/boxplot
    plotBarBoxMetrics(hAxB,xi,p,pStr,pP,yL);      
    
    % updates the axis properties
    set(hAxB,'xticklabel',[],'xlim',xLim,'linewidth',1.5,'UserData',1);

    % formats the single axis
    formatPlotAxis(hAxB,pF,1);                     
    
    % turns the grid on (if specified)
    if (pP.plotGrid); set(hAxB(1),'ygrid','on'); end      
else
    % removes the label indices
    [pF.xLabel.ind,pF.yLabel.ind] = deal(NaN);
    
    % loops through all the indices plotting the solutions
    for j = 1:nApp
        % initialises the subplot axes
        [i,hAxB{j}] = deal(ind(j),createSubPlotAxes(hP,[m,n],j));
                
        % plots the background image
        if pP.pltDN
            fill(xFillD(ix),yFill(iy),'y','FaceAlpha',fAlpha',...
                'tag','hDN','HitTest','off')
            fill(xFillN(ix),yFill(iy),'k','FaceAlpha',fAlpha',...
                'tag','hDN','HitTest','off')
        end      

        % creates the bar graph/boxplot
        if (strcmp(pStr,'Tw'))
            plotBarBoxMetrics(hAxB{j},xi,p(j),pStr,pP,yL);
            if (j == 1); yLim = 60; end
        else
            yLim = max(yLim,plotBarBoxMetrics(hAxB{j},xi,p(j),pStr,pP,yL));    
        end        

        % updates the axis properties
        set(hAxB{j},'xticklabel',[],'linewidth',0.5);
        formatPlotAxis(hAxB{j},pF,i);                
        
        % turns the grid on (if specified)
        if (pP.plotGrid); set(hAxB{j},'ygrid','on'); end         
    end    
end
    
% ensures all the y-axis for the sleep bouts are the same
if (nGrp > 1)
    if (~isempty(hAxB{1}))
        for i = 1:length(hAxB)
            if (isempty(yL))
                yLimNw = setStandardYAxis(hAxB{i}(1),[],nTick,yLim);
            else
                yLimNw = yL(2);
            end
        end

        % resets the height of the day/night background to the overall limits
        if (pP.pltDN)
            for i = 1:length(hAxB{1})
                hFill = cellfun(@(x)(findobj(x(i),'tag','hDN')),hAxB,'un',0);
                yDataNw = [0 yLimNw*[1 1] 0];
                cellfun(@(x)(set(x,'yData',yDataNw)),hFill);
            end
        end
    else
        % retrieves the fill handles
        hFill = findall(get(hAx,'parent'),'tag','hDN');
        set(hFill,'yData',[0 60*[1 1] 0]);    
    end
end

% sets the non-aligned x/y labels
formatMultiXYLabels(hAxB,pF,[m,n]);    

% resets the axis positions
if (nGrp == 1)
    % case is for a single group
    resetAxesPos(hAxB,m,n,[0,0.04]);
    resetObjPos(hAxB,'Bottom',0.02,1);
    
    % sets the group strings
    setGroupString(hAxB,pF,xi,snTot(1).iMov.pInfo.gName(ind),tRot,-0.005);
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