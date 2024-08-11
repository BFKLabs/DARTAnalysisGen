% --- plots the daily average speed of the flies over the experiment over a 
%     specified time bin duration (long experiment only)
function pData = AvgSpeedPlotDN(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Population Movement (Long)';
pData.Type = {'Pop','Multi'};
pData.fType = [1 1 3 1];
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
    pData.canComb = length(snTot(1).iMov.ok) > 1;
    pData.hasSP = true;
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
rI.Dur = 'Long';
rI.Shape = 'None';
rI.Stim = 'None';
rI.Spec = 'None';

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 1;
cP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'Number',60,'tBin','Time Bin Size (sec)',[5 3600 false]);

% sets the tool-tip strings
cP(1).TTstr = 'Duration over which the population speed is averaged';

% --- initialises the plotting parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 5;
pP = setParaFields(nPara);

% sets the parameter list for 
pList = {'SEM','Min/Max'};

% sets the tab list names
a = {'1 - General','2 - Signal Error'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},'Number',1,'lWid','Plot Line Width',[0.1 10 false]);
pP(2) = setParaFields(a{1},'Boolean',0,'isZeitG','Use Zeitgeiber Time Format');
pP(3) = setParaFields(a{1},'Boolean',0,'plotGrid','Show Plot Axis Gridlines');
pP(4) = setParaFields(a{2},'Boolean',0,'pltErr','Plot Signal Error');
pP(5) = setParaFields(a{2},'List',{1,pList},'errType','Signal Error Type',[],{4,2});

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
    pF.Title(i).String = snTot(1).iMov.pInfo.gName{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);

% sets the independent variable fields
oP = addXVarField(oP,'Time','T','Time');

% appends the output metric parameter fields
oP = addYVarField(oP,'Speed','V',[],5,{'T'},1);
oP = addYVarField(oP,'Speed (Mean)','V_mn',[],4,{'T'},1);
oP = addYVarField(oP,'Speed (SEM)','V_sem',[],4,{'T'},1);
oP = addYVarField(oP,'Speed (Min)','V_min',[],4,{'T'},1);
oP = addYVarField(oP,'Speed (Max)','V_max',[],4,{'T'},1);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% global variables
global tDay

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
sP = retParaStruct(dcObj.pData.sP);

% sets the common class fields
dcObj.yName = 'Average Speed';
dcObj.xName = 'Time';
dcObj.pType = 'Trace';
dcObj.tUnits = 'Hours';
dcObj.yUnits = 'mm/sec';
dcObj.combFig = sP.Sub.isComb;
dcObj.tDay0 = [0,0,0,tDay,0,0];
dcObj.grpName = dcObj.pData.appName(sP.Sub.isPlot);

% sets up the data cursor string
dTxt = dcObj.setupCursorString();

% ----------------------------------------------------------------------- %
% ---                       CALCULATION FUNCTION                      --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [plotD,ok] = calcFunc(snTot,pData,gPara,cP)

% initialises the calculation parameters (if not already initialised)
if (nargin == 3)
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
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);

% sets the binned time vector
T = ((cP.tBin/2):(cP.tBin):convertTime(1,'day','sec'))';

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,'T',T,'V',[],'V_mn',[],...
                                 'V_sem',[],'V_min',[],'V_max',[],...
                                 'tBin',cP.tBin,'movType',cP.movType,...
                                 'hDay',gPara.TdayC);

% creates the waitbar figure
wStr = {'Overall Progress','Average Speed Calculations'};
wOfs = (nExp > 1); wStr = wStr((2-wOfs):end);
h = ProgBar(wStr,'Average Speed Calculations');

% determines the binned indices
h.Update(1,'Determining Time Bin Indices...',1/(2+nApp));

% memory allocation
for j = 1:nApp
    plotD(j).V = cellfun(@(x)(NaN(length(T),1)),plotD(j).V,'un',0);
end

% ----------------------------- %
% --- VELOCITY CALCULATIONS --- %
% ----------------------------- %

% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure (if more than one solution file)
    if wOfs > 0
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
        if h.Update(1,wStrNw,i/(1+nExp))
            [plotD,ok] = deal([],false);
            return
        end
    end
    
    % determines the index offset
    Ttot = cell2mat(snTot(i).T);
    flyok = snTot(i).iMov.flyok;    
    
    % determines the binned indices (for time length, tBin) and determines the
    % bins which has at least two time points
    h.Update(1+wOfs,'Determining Time Bin Indices...',0.25);
    indB = detTimeBinIndices(Ttot,cP.tBin,1);
    
    % determines the time group index array
    [Tnw,indB] = setBinnedTimeArray(Ttot,indB,cP.tBin);
    indD = detTimeGroupIndexArray(Tnw,snTot(i).iExpt(1),cP.tBin,cP.Tgrp0);
        
    % calculates the binned activity 
    for j = 1:nApp
        % updates the waitbar figure
        wStrNw = sprintf(['Calculating Binned Activity (Region ',...
                          '%i of %i)'],j,nApp);
        if h.Update(1+wOfs,wStrNw,0.50*(1+(j/(nApp+1))))
            % if the user cancelled, then exit the function
            [plotD,ok] = deal([],false);
            return            
        end     
    
        % only calculate if values exist...
        if ~isempty(snTot(i).Px{j})
            % calculates the binned fly movement speed/midline crossings
            V = calcBinnedFlyMovement(snTot(i),Ttot,indB,cP,j,flyok{j});    
        else
            % if no valid flies then set a NaN array
            V = repmat({NaN(1,length(flyok{j}))},length(indB),1);                     
        end        

        % sets the raw velocity values into the plotting data struct          
        Vr = num2cell(cell2mat(V),1);
        if ~isempty(Vr)
            for iD = 1:size(indD,2)
                % sets the velocity values (for all flies) for each day
                iiN = ~isnan(indD(:,iD));            
                for iFly = find(snTot(i).iMov.flyok{j}(:)')
                    plotD(j).V{iD,iFly,i}(iiN) = Vr{iFly}(indD(iiN,iD));
                end
            end
        end
    end
end                    
    
% updates the waitbar figure
for j = 1:nApp
    wStrNw = sprintf(['Calculating Population Movement (Region ',...
                      '%i of %i)'],j,nApp);
    if h.Update(1+wOfs,wStrNw,(j+1)/(2+nApp))
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    end                 
                          
    % calculates the metric statistics (removes any infinite SEM values)
    plotD(j) = calcMetricStats(plotD(j),'V');
    plotD(j).V_sem(isinf(plotD(j).V_sem)) = NaN;     
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
cP = retParaStruct(pData.cP);
sP = retParaStruct(pData.sP);
pF = pData.pF;

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the plotting indices and subplot indices
[ind,m,n] = deal(find(sP.Sub.isPlot),sP.Sub.nRow,sP.Sub.nCol);
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% sets the time multiplier (from seconds to hours)
sP.xLim = [0 convertTime(24,'hrs','sec')];
[movType,tBin] = deal(p(1).movType,p(1).tBin);

% sets the alignment flags and time multipliers
[isAlign,absTime] = deal(false,true);
tMlt = convertTime(1,'sec','hrs');   

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the subplot titles
if (sP.Sub.isComb)
    % case is combining, so remove the titles
    [pF.xLabel.ind,pF.yLabel.ind] = deal(1); 
    [pF.Legend.String,m,n] = deal(snTot(1).iMov.pInfo.gName(ind),1,1);
    for i = 1:length(pF.Title); pF.Title(i).String = ''; end
else
    % resets the labels
    [pF.xLabel.ind,pF.yLabel.ind] = deal(NaN); 
    [pF.xLabel.String,pF.yLabel.String] = deal(''); 
    
    % sets the axis titles
    for i = 1:length(pF.Title)
        pF.Title(i).String = snTot(1).iMov.pInfo.gName{i};         
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
    % case is using the absolute time axis
    if (pP.isZeitG)            
        % case is using Zeitgeiber time
        pF.xLabel.String = 'Zeitgeiber Time';
    else
        % if plotting day/night, then set absolute time
        pF.xLabel.String = 'Time'; 
    end
end

% determines if the signals can/are aligned
if isfield(cP,'isAlign')
    if cP.isAlign
        % if the signals are aligned, then remove the other flags
        isAlign = true;
        [pP.pltDN,pP.isZeitG,absTime] = deal(false);   
        
        % determines if experiment is very short (< 2 hours). if so, then
        % user minutes as a time scale (instead of hours)
        Tf = ceil(cellfun(@(x)(x{end}(end)),field2cell(snTot,'T')));
        if max(convertTime(Tf,'sec','hours') < 2)
            tMlt = convertTime(1,'sec','mins');     
            pF.xLabel.String = 'Time (min)'; 
        end
    end    
end

% sets the output data
pData.pF = pF;

% retrieves the formatting struct
if isempty(m); szMx = 1; else; szMx = max([m n]); end
pF = retFormatStruct(pF,szMx);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% memory allocation
hPlot = cell(nApp,1);
if pP.pltErr
    % calculates the overall limit (mean + error)
    switch pP.errType
        case ('SEM') % error type is SEM
            yLim = detOverallLimit(cellfun(@(x)(max(cell2mat(...
                field2cell(x,'V_mn'))+cell2mat(field2cell(x,'V_sem')))),plotD));
        case ('Min/Max') % error type is min/max
            yLim = detOverallLimit(cellfun(@(x)(max(cell2mat(...
                field2cell(x,'V_max')))),plotD));
    end
    
    % if not SEM, then calculate maximum from the mean
    if isnan(yLim)
        yLim = max(cellfun(@(x)(detOverallLimit(field2cell(x,'V_mn'))),plotD));    
    end
else
    % calculates the overall limit (mean)
    yLim = max(cellfun(@(x)(detOverallLimit(field2cell(x,'V_mn'))),plotD));
end
    
% sets the trace colours
[col,wStr] = deal(num2cell(distinguishable_colors(nApp,'w'),2),[]);
hAx = cell(1+(nApp-1)*(~sP.Sub.isComb),1);    

% loops through all the subplots 
for j = 1:nApp
    % sets the actual plot index
    i = ind(j);

    % sets the subplot to be plotted on and updates the properties
    if sP.Sub.isComb
        % case is combining, so use the main axis
        colNw = col{j};
        if (j == 1)
            [hAx{j},hAxNw] = deal(createSubPlotAxes(hP)); 
            hold(hAx{j},'on'); axis(hAxNw,'on')
        end

        % adds the day/night (if not aligned)
        if ~isAlign && (j == 1); plotSingleDNGraph(yLim,p(j)); end        
    else
        % otherwise, create a seperate subplot
        colNw = col{1};
        [hAx{j},hAxNw] = deal(createSubPlotAxes(hP,[m,n],j)); 
        hold(hAx{j},'on');
        
        % adds the day/night (if not aligned)
        if ~isAlign; plotSingleDNGraph(yLim,p(j)); end           
    end
    
    % turns the axis box on
    set(hAxNw,'linewidth',0.5,'box','on')

    % if there is an overlap, then plot the SEM traces
    Tplt = (p(j).T-p(j).T(1))*tMlt + convertTime(tBin/2,'sec','hrs');
    if pP.pltErr
        switch pP.errType
            case ('SEM')
                % sets the error signal and y-axis limit                    
                Yerr = p(j).V_sem;
            case ('Min/Max')
                % sets the error signal and y-axis limit                    
                Yerr = [p(j).V_min p(j).V_max];
        end

        % plots the SEM signal
        plotSignalSEM(p(j).V_mn,Yerr,Tplt,colNw,0.5)             
    end

    % plots the traces
    hPlotNw = plot(hAxNw,Tplt,p(j).V_mn,'color',colNw,...
        'linewidth',pP.lWid,'UserData',j);
    if sP.Sub.isComb; hPlot{j} = hPlotNw; end   

    % turns the grid on/off depending on the flag values
    if pP.plotGrid
        % grid is turned on
        grid on
    else
        % grid is turned off
        grid off
    end           
    
    % sets the x/y axis limits         
    if isAlign
        Tf = ceil(cellfun(@(x)(x{end}(end)),field2cell(snTot,'T')));
        set(hAxNw,'xlim',[0 max(Tf)*tMlt]+0.05*[-1 1]);
    else
        if ~isnan(yLim)
            set(hAxNw,'xlim',sP.xLim*tMlt,'ylim',[0 yLim])        
        else
            set(hAxNw,'xlim',sP.xLim*tMlt)                        
        end
    end
    
    % case is using the absolute time axis
    Tax = convertTime([0 24],'hrs','sec');
    if pP.isZeitG
        % case is using Zeitgeiber time
        setZeitGTimeAxis(hAxNw,Tax);        
    else
        % if plotting day/night, then set absolute time
        if absTime
            setAbsTimeAxis(hAxNw,Tax);  
        end
    end       
    
    % formats the plot axis labels
    expandAxesLim(hAxNw);
    formatPlotAxis(hAxNw,pF,i);
end

% --- PLOT AXES REFORMATTING --- %
% ------------------------------ %

% sets the non-aligned x/y labels
if ~sP.Sub.isComb
    formatMultiXYLabels(hAx,pF,[m,n]);
end

% formats and resets the axis positions
resetAxesPos(hAx,m,n); 

% updates the figure if plotting a combined figure
if sP.Sub.isComb
    % creates the legend object
    hLg = createLegendObj(hPlot,pF.Legend);

    % resets the legend position           
    pLg = get(hLg,'position');
    set(hLg,'Position',[(1-pLg(3)),0.5*(1-pLg(4)),pLg(3:4)])   
    resetLegendPos(hLg,hAx)      
    
    % resets the axes position        
    set(hAxNw,'ylim',get(hAxNw,'ylim'))
else
    % resets the axis maximum limits
    yLimMx = max(cellfun(@(x)(max(get(x,'ylim'))),hAx));
    cellfun(@(x)(set(x,'ylim',[0 yLimMx])),hAx)        
end

% if there is a warning string, then show it to screen
if ~isempty(wStr)
    waitfor(warndlg(wStr))
end

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)
