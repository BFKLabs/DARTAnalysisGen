% --- plots the average speed of the flies over the duration of an 
%     experiment for a specified time bin duration
function pData = AvgSpeedPlot(snTot)

% initialise the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);

% sets the function name/type
pData.Name = 'Population Movement (Full Experiment)';
pData.Type = {'Pop'};
pData.fType = [1 1 1 1];
pData.rI = initFuncReqInfo(pData);
pData.dcFunc = @dataCursorFunc;

% initialises the other fields  (if input argument provided)
if nargin == 1   
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
    pData.canComb = length(snTot(1).iMov.ok) > 1;
    [pData.hasSP,pData.hasTime] = deal(true);
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
nPara = 2;
cP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'Boolean',0,...
            'calcInst','Calculate Instantaneous Speed');
cP(2) = setParaFields(a{1},'Number',60,'tBin',...
            'Activity Calculation Duration (sec)',[1 3600 false],{1,1});

% sets the tool-tip strings
cP(1).TTstr = 'Calculates instantaneous speed intead of binned speed';
cP(2).TTstr = 'Duration over which the population speed is averaged';

% adds the unique motor parameters
cP = addUniqueMotorPara(cP,snTot);

% --- initialises the plotting parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
isLong = ~detIfShortExpt(field2cell(snTot,'T'));
Ts = getMotorFiringTimes(snTot.stimP);
hasStim = ~isempty(Ts);
nPara = 2 + hasStim + 2*isLong;
pP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General','2 - Signal Error'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},...
                    'Number',1,'lWid','Plot Line Width',[0.1 10 false]);
pP(2) = setParaFields(a{2},'Boolean',0,'pltErr','Plot Signal Error');

% adds the show stimuli marker (if there were stimuli)
if hasStim
    pP(3) = setParaFields(a{1},...
                'Boolean',1,'showStim','Show Stimuli Markers');
end

% sets the parameters for the day/night and time format
if isLong
    pP(3+hasStim) = setParaFields(a{1},...
                'Boolean',1,'pltDN','Plot Day/Night Background Image');
    pP(4+hasStim) = setParaFields(a{1},...
                'Boolean',0,'isZeitG','Use Zeitgeiber Time Format');
end
    
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

% sets the independent variable fields
oP = addXVarField(oP,'Time','T','Time');

% sets the dependent variable fields
oP = addYVarField(oP,'Speed','V',[],5,{'T'},1);
oP = addYVarField(oP,'Speed (Mean)','V_mn',[],4,{'T'},1);
oP = addYVarField(oP,'Speed (SEM)','V_sem',[],4,{'T'},1);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
sP = retParaStruct(dcObj.pData.sP);

% sets the common class fields
dcObj.yName = 'Average Speed';
dcObj.xName = 'Time';
dcObj.pType = 'Trace';
dcObj.yUnits = 'mm/sec';
dcObj.combFig = sP.Sub.isComb;
dcObj.grpName = dcObj.pData.appName(sP.Sub.isPlot);

% sets up the time scale values
Tend = dcObj.plotD{1}(1).T(end);
[~,dcObj.tUnits] = getTimeScale(Tend);

% sets up the data cursor string
dTxt = dcObj.setupCursorString();

% ----------------------------------------------------------------------- %
% ---                       CALCULATION FUNCTION                      --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [plotD,ok] = calcFunc(snTot,pData,gPara,cP)

% initialises the calculation parameters (if not already initialised)
if nargin == 3
    % retrieves the parameter struct
    cP = retParaStruct(pData.cP,gPara);

    % sets the movement type (based on the global parameters)
    if strcmp(gPara.movType,'Absolute Location')
        cP.movType = 'Absolute Speed';
    end
end
    
% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
[nApp,ok] = deal(length(snTot.iMov.ok),true);
Ttot = cell2mat(snTot.T);
flyok = snTot.iMov.flyok;

% sets up the time vector and time binning arrays (if required)
if cP.calcInst
    % case is calculating instantaeous speed
    T = Ttot(2:end);
    [dT,isOK] = deal(diff(Ttot),true(size(T)));
    [TD,indD] = deal(convertTime(T/60,'min','sec'),(1:length(T))');
    
else
    % determines the binned indices (for time length, tBin) and determines 
    % the bins which has at least two time points
    indB = detTimeBinIndices(Ttot,cP.tBin);
    [T,indB] = setBinnedTimeArray(Ttot,indB,cP.tBin);

    % determines the 
    indD = detTimeGroupIndexArray(T,snTot.iExpt(1),cP.tBin,cP.Tgrp0);
    isOK = find(~isnan(indD));
    TD = convertTime((cP.tBin*(1:length(T))' - cP.tBin/2)/60,'min','sec');
end

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,'movType',cP.movType,...
                            'Tf',[],'T',TD,'V',[],'V_mn',[],'V_sem',[],...
                            'TInfo',snTot.iExpt.Timing);

% creates the waitbar figure
wStr = {'Overall Progress'};
h = ProgBar(wStr,'Average Speed Calculations');

% -------------------------- %
% --- SPEED CALCULATIONS --- %
% -------------------------- %

% determines the binned indices
h.Update(1,'Determining Time Bin Indices...',1/(2+nApp));

% calculates the movement traces for all the selected apparatus
for j = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf(['Calculating Population Movement (Region ',...
                      '%i of %i)'],j,nApp);
    if h.Update(1,wStrNw,(j+1)/(2+nApp))
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    end
    
    % determines if there any valid flies in the group
    if any(flyok{j})
        % only calculate if values exist...    
        if ~isempty(snTot.Px{j})
            if cP.calcInst
                % calculates the instantaneous speed
                V = calcInstantSpeed(snTot,dT,j);
            
            else
                % calculates binned fly movement speed/midline crossings
                V = calcBinnedFlyMovement(snTot,Ttot,indB,cP,j,flyok{j});
            end
            
        else
            % otherwise, set a NaN array
            V = repmat({NaN(1,length(flyok{j}))},length(indB),1); 
        end
    else
        % if no valid flies then set a NaN array
        V = repmat({NaN(1,length(flyok{j}))},length(indB),1); 
    end
    
    % removes any empty speed regions from the time array
    if j == 1; isNE = ~cellfun('isempty',V); end
    plotD(j).Tf = T(isNE);
    
    % sets the raw velocity values            
    plotD(j) = setRawDataValues(plotD(j),snTot,V,indD,'V',1,j,1);            

    % calculates the statistic metrics
    plotD(j) = calcMetricStats(plotD(j),'V',4,1);
    if length(plotD(j).T) ~= length(isOK)
        plotD(j).T = plotD(j).T(indD(isOK));
        plotD(j).Tf = plotD(j).Tf(indD(isOK));
    end
end

% closes the waitbar figure
if ~h.Update(1,'Speed Calculations Complete!',1)
    h.closeProgBar()
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

% sets any missing fields
if ~isfield(pP,'showStim'); pP.showStim = false; end

% retrieves the other calculation parameters (if they exist)
[devType,chType] = deal([]);
if isfield(cP,'devType'); devType = cP.devType; end
if isfield(cP,'chType'); chType = cP.chType; end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the plotting indices and subplot indices
[ind,m,n] = deal(find(sP.Sub.isPlot),sP.Sub.nRow,sP.Sub.nCol);
nApp = length(ind); if nApp == 0; return; end
p = plotD{1}(ind);

% if the DN parameters are not, then set default values
if ~isfield(pP,'pltDN')
    [pP.pltDN,pP.isZeitG,absTime] = deal(false);    
else
    absTime = true;
end

% sets the movement type and time multiplier
movType = p(1).movType;
tMlt = getTimeScale(snTot.T{end}(end));

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the subplot titles
if sP.Sub.isComb
    % case is combining, so remove the titles    
    [pF.xLabel.ind,pF.yLabel.ind] = deal(1);  
    [pF.Legend.String,m,n] = deal(snTot.iMov.pInfo.gName(ind),1,1);
    for i = 1:length(pF.Title); pF.Title(i).String = ''; end    
else    
    % resets the labels
    [pF.xLabel.ind,pF.yLabel.ind] = deal(NaN);  
    [pF.xLabel.String,pF.yLabel.String] = deal(''); 
    
    % sets the axis titles
    for i = 1:length(pF.Title)
        pF.Title(i).String = snTot.iMov.pInfo.gName{i};         
    end
end

% resets the y-label string (if using mid-line crossings)
if isempty(pF.yLabel.String)
    if strcmp(movType,'Midline Crossing')
        pF.yLabel.String = sprintf('Beam Crosses (count min^{-1})');
    else
        pF.yLabel.String = sprintf('Speed (mm sec^{-1})');
    end
end
    
% sets the time axis properties
if isempty(pF.xLabel.String)
    if absTime
        % case is using the absolute time axis
        if pP.isZeitG
            % case is using Zeitgeiber time
            pF.xLabel.String = 'Zeitgeiber Time';
        else
            % case is not using Zeitgeiber time
            if pP.pltDN
                % if plotting day/night, then set absolute time
                pF.xLabel.String = 'Time'; 
            else
                % otherwise, reset the xlabel string
                pF.xLabel.String = 'Time (hours)'; 
            end
        end
    else
        % otherwise, reset the xlabel string
        pF.xLabel.String = 'Time (mins)';
    end            
end

% sets the output data
pData.pF = pF;

% retrieves the formatting struct
if isempty(m); szMx = 1; else; szMx = max([m n]); end
pF = retFormatStruct(pData.pF,szMx);

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% memory allocation
hPlot = cell(nApp,1);
if pP.pltErr
    % calculates the overall limit (mean + SEM)
    yLim = detOverallLimit(cellfun(@(x)(max(cell2mat(...
        field2cell(x,'V_mn'))+cell2mat(field2cell(x,'V_sem')))),plotD));
    
    % if not SEM, then calculate maximum from the mean
    if isnan(yLim)
        yLim = max(cellfun(@(x)(detOverallLimit(field2cell(x,'V_mn'))),plotD));    
    end
else
    % calculates the overall limit (mean)
    yLim = max(cellfun(@(x)(detOverallLimit(field2cell(x,'V_mn'))),plotD));
end

% sets the stimulus markers
if pP.showStim
    Ts = getMotorFiringTimes(snTot.stimP,devType,chType);
    if iscell(Ts); Ts = Ts{1}; end
    [Tstim,Ystim] = deal(repmat(Ts,1,2)',repmat([0 yLim],length(Ts),1)');
end
    
% plots the day/night bands (if required)
yOfs = 0;
if pP.pltDN
    if sP.Sub.isComb
        hAx = {plotDayNightGraph(hP,snTot,p(1).Tf*tMlt,ceil(yLim+yOfs),1,[m,n],tMlt)};
    else
        hAx = plotDayNightGraph(hP,snTot,p(1).Tf*tMlt,ceil(yLim+yOfs),ind,[m,n],tMlt);
        if (m*n) == 1; hAx = {hAx}; end
    end
else
    if sP.Sub.isComb
        hAx = {createSubPlotAxes(hP)};
    else
        hAx = arrayfun(@(x)(createSubPlotAxes(hP,[m,n],x)),1:nApp,'un',0);
    end
end 

% sets the trace colours
col = num2cell(distinguishable_colors(nApp,'w'),2);

% loops through all the subplots 
for j = 1:nApp
    % sets the actual plot index
    i = ind(j);
        
    % sets the subplot to be plotted on and updates the properties
    if sP.Sub.isComb
        % case is combining, so use the main axis
        [colNw,k] = deal(col{j},1);
        if j == 1; axis(hAx{k},'on'); end
    else
        % otherwise, create a seperate subplot
        [colNw,k] = deal(col{1},j);
    end        
    
    % turns the axis box on  
    hold(hAx{k},'on');
    set(hAx{k},'linewidth',1.5,'box','on')

    % plots the SEM error signal (if required)
    if pP.pltErr 
        set(gcf,'CurrentAxes',hAx{k})
        plotSignalSEM(p(j).V_mn+yOfs,p(j).V_sem,p(j).Tf*tMlt,colNw,0.5)     
    end        
        
    % plots the stimulus markers
    if (i == 1) || ~sP.Sub.isComb
        if pP.showStim && exist('Tstim','var')
            plot(hAx{k},Tstim*tMlt,Ystim+yOfs,'k:','linewidth',...
                pP.lWid,'tag','hStim','HitTest','off'); 
        end    

        % sets the x/y axis limits
        if range(sP.xLim) == 0
            set(hAx{k},'ylim',[0 yLim]+yOfs,'box','off')            
        else
            ylimT = [0 yLim]+yOfs;
            set(hAx{k},'xlim',sP.xLim*tMlt,'ylim',ylimT,'box','off')            
        end     

        % sets the time axis properties
        if pP.isZeitG    
            % case is using Zeitgeiber time
            setZeitGTimeAxis(hAx{k},p(j).Tf,snTot);        
        else
            % if plotting day/night, then set absolute time
            setAbsTimeAxis(hAx{k},p(j).Tf,snTot);        
        end        
        
        % formats the plot axis
        formatPlotAxis(hAx{k},pF,i); 
        set(hAx{k},'UserData',k);
        axis(hAx{k},'on')
    end
    
    % plots the traces (time scaled to hours)  
    hPlotNw = plotFullSignal(hAx{k},p(j).Tf*tMlt,p(j).V_mn+yOfs);
    set(hPlotNw,'color',colNw,'linewidth',pP.lWid,'UserData',j);    
    if sP.Sub.isComb; hPlot{j} = hPlotNw; end              
end

% ------------------------------ %
% --- PLOT AXES REFORMATTING --- %
% ------------------------------ %

% sets the non-aligned x/y labels
if ~sP.Sub.isComb
    formatMultiXYLabels(hAx,pF,[m,n]);
end

% updates the figure if plotting a combined figure
if sP.Sub.isComb
    % creates the legend object
    [hLg,yLimMx] = deal(createLegendObj(hPlot,pF.Legend),yLim);           
else
    % resets the axis maximum limits
    yLimMx = max(cellfun(@(x)(max(get(x,'ylim'))),hAx));
    cellfun(@(x)(set(x,'ylim',[0+yOfs yLimMx])),hAx)         
end

% sets up the y-axis properties (for each axes)
yLimMx = cellfun(@(x)(setStandardYAxis(x,[],5,yLimMx,0)),hAx);
for i = 1:length(hAx)
    % resets the stimuli markers
    if pP.showStim
        set(findall(hAx{i},'tag','hStim'),'yData',[0 yLimMx(1)]);
    end
    
    % resets the day/night background (if plotted)
    if pP.pltDN
        hDN = findall(hAx{i},'tag','hDN');
        yData = get(hDN(1),'yData');
        yData(yData>mean(yData)) = yLimMx(i);
        set(hDN,'yData',yData)
    end
end

% formats and resets the axis positions
resetAxesPos(hAx,m,n); 

% resets the legend location
if sP.Sub.isComb
    resetLegendAxisPos(hAx{1},hLg,[1 1],-0.005)
end

% ----------------------------------------------------------------------- %
% ---                     MISCELLANEOUS FUNCTION                      --- %
% ----------------------------------------------------------------------- %

% --- calculates the instantaneous speed
function V = calcInstantSpeed(snTot,dT,ind)

% sets calculates the squared displacement
dP2 = diff(snTot.Px{ind},[],1).^2;
if ~isempty(snTot.Py)
    % case is a 2D experiments
    dP2 = dP2 + diff(snTot.Py{ind},[],1).^2; 
end

% calculates the instantaneous speed
V = num2cell(sqrt(dP2)./dT,2);

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)