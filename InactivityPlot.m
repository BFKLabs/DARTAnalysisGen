% --- plots the daily fly proportional activity/sleep mins per hour 
%    (full experiment)
function pData = InactivityPlot(snTot)

% initialise the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);

% sets the function name/type
pData.Name = 'Temporal Inactivity (Full Experiment)';
pData.Type = {'Pop','Multi'};
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
nPara = 1;
cP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'Number',...
    10,'tMove','Inactivity Detection Duration (sec)',[1 20 false]);

% sets the tool-tip strings
cP(1).TTstr = 'Time period over which inactivity proportion is calculated';

% --- initialises the plotting parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
isLong = ~detIfShortExpt(field2cell(snTot,'T'));
nPara = 7 + isLong;
pP = setParaFields(nPara);

% sets the parameter lists
pList = {'Hourly Sleep Duration','Inactivity Proportion'};

% sets the tab list names
a = {'1 - General','2 - Signal Error','3 - Smoothing'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},...
    'List',{1,pList},'inactType','Inactivity Metric');    
pP(2) = setParaFields(a{1},...
    'Number',1,'lWid','Plot Line Width',[0.1 10 false]);
pP(3) = setParaFields(a{1},...
    'Boolean',0,'isZeitG','Use Zeitgeiber Time Format');
pP(4) = setParaFields(a{1},...
    'Boolean',0,'plotGrid','Show Plot Axis Gridlines');
pP(5) = setParaFields(a{2},...
    'Boolean',0,'plotErr','Plot Signal Error');
pP(6) = setParaFields(a{3},...
    'Boolean',0,'useSm','Use Signal Smoothing');
pP(7) = setParaFields(a{3},...
    'Number',11,'nSm','Smoothing Window Size',[3 50 true],{6,2});

% sets the parameters for the day/night and time format
if isLong
    pP(8) = setParaFields(a{1},...
                'Boolean',1,'pltDN','Plot Day/Night Background Image');
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

% appends the output metric parameter fields
oP = addYVarField(oP,'Inactivity','I',[],5,{'T'},1);
oP = addYVarField(oP,'Inactivity (Mean)','I_mn',[],4,{'T'},1);
oP = addYVarField(oP,'Inactivity (SEM)','I_sem',[],4,{'T'},1);
oP = addYVarField(oP,'Inactivity (Min)','I_min',[],4,{'T'},1);
oP = addYVarField(oP,'Inactivity (Max)','I_max',[],4,{'T'},1);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% global variables
global tDay

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
cP = retParaStruct(dcObj.pData.cP);
sP = retParaStruct(dcObj.pData.sP);

% sets the common class fields
dcObj.xName = 'Time';
dcObj.pType = 'Trace';
dcObj.tUnits = 'Hours';
dcObj.combFig = sP.Sub.isComb;
dcObj.tDay0 = [0,0,0,tDay,0,0];
dcObj.grpName = dcObj.pData.appName(sP.Sub.isPlot);
dcObj.yUnits = '%';        

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
        cP.movType = 'Absolute Range';
    end
end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
nExp = length(snTot);
nApp = length(snTot(1).iMov.ok);
[snTot,Ttot] = setupTimeVectors(snTot);

% creates the waitbar figure
wOfs = (nExp > 1); 
wStr0 = {'Setting Up Time Bin Arrays...','Inactivity Calculations'};
wStr = wStr0(1:(1+wOfs));
h = ProgBar(wStr,'Inactivity Proportion Calculations');
pause(0.05);

% sets up the time bin data
[T,indB,indD,useInt,ok] = setupTimeBinData(snTot,Ttot,cP,h);
if ~ok
    % if the user cancelled, then exit
    plotD = [];
    return
end

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,'T',T,'I',[],'I_mn',[],...
                             'I_sem',[],'I_min',[],'I_max',[],...                             
                             'hDay',gPara.TdayC); 

% other initialisations
nT = length(T);

% ------------------------------- %
% --- INACTIVITY CALCULATIONS --- %
% ------------------------------- %
                         
% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure (if more than one solution file)
    if wOfs > 0
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
        if h.Update(1,wStrNw,i/(1+nExp))
            [plotD,ok] = deal([],false);
            return
        else
            h.Update(2,wStr{2},0);
            pause(0.05);
        end
    end

    % determines the index offset
    fOK = snTot(i).iMov.flyok;
    
    % sets the use interpolation flag for the experiment
    if useInt && (nT ~= (length(Ttot{i}) - 1))
        snTot(i) = interpFlyCoords(snTot(i),[0;T],Ttot{i});
    end    
    
    % calculates the binned activity 
    for j = 1:nApp
        % if there are no coordinates, then continue
        if isempty(snTot(i).Px{j})
            continue
        end
        
        % calculates the inactivity metrics based on the calculation type        
        % updates the waitbar figure
        wStrNw = sprintf(['Calculating Proportional Inactivity ',...
                          '(Region %i of %i)'],j,nApp);
        if h.Update(1+wOfs,wStrNw,0.50*(1+(j/(nApp+1))))
            % if the user cancelled, then exit the function
            [plotD,ok] = deal([],false);
            return                       
        end

        % calculates the binned activity
        Anw = calcBinnedActivity(snTot(i),Ttot{i},indB{i},cP,j,fOK{j});
        Inw = cellfun(@(x)(~x),Anw,'un',0);

        % sets the raw values into the storage array
        plotD(j) = setRawDataValues(...
                    plotD(j),snTot(i),Inw,indD{i},'I',i,j,1);           
    end
end

% calculates the final metrics
for j = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf(['Calculating Proportional Inactivity (Region ',... 
                      '%i of %i)'],j,nApp);
    if h.Update(1+wOfs,wStrNw,(j+1)/(2+nApp))
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    end                         
                                           
    % calculates the statistical metrics
    plotD(j) = calcMetricStats(plotD(j),'I',4,1);
end

% closes the waitbar figure
if ~h.Update(1,'Inactivity Calculations Complete!',1)
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

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% sets the other plot properties
tMlt = getTimeScale(p(1).T);
isInact = strcmp(pP.inactType,'Inactivity Proportion');

% if the DN parameters are not, then set default values
if ~isfield(pP,'pltDN')
    pP.pltDN = false;
end

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %
    
% sets the subplot titles
if sP.Sub.isComb
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

% sets the x/y labels
pF.yLabel.String = '% Inactivity'; 

% sets the time axis properties
if isempty(pF.xLabel.String)
    % case is using the absolute time axis
    if pP.isZeitG          
        % case is using Zeitgeiber time
        pF.xLabel.String = 'Zeitgeiber Time';
    else
        % if plotting day/night, then set absolute time
        pF.xLabel.String = 'Time'; 
    end
end

% sets the time axis properties (based on metric type)
if isInact
    % case is inactivity proportion
    [yLim,pY] = deal(100);
    
else
    % case is the hourly sleep duration
    [yLim,pY] = deal(60);
    pF.yLabel.String = sprintf('Sleep (Min/Hour)'); 
end

% retrieves the formatting struct
if isempty(m); szMx = 1; else; szMx = max([m n]); end
pF = retFormatStruct(pF,szMx);
    
% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% memory allocation
hPlot = cell(nApp,1);   

% plots the day/night bands (if required)
if pP.pltDN
    if sP.Sub.isComb
        hAx = {plotDayNightGraph(hP,...
            snTot(1),p(1).T*tMlt,yLim,1,[m,n],tMlt)};
    else
        hAx = plotDayNightGraph(hP,...
            snTot(1),p(1).T*tMlt,yLim,ind,[m,n],tMlt);
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
    Tplt = p(j).T*tMlt;    
    
    % sets the subplot to be plotted on and updates the properties
    if sP.Sub.isComb
        % case is combining, so use the main axis
        [colNw,k] = deal(col{j},1);
    else
        % otherwise, create a seperate subplot
        [colNw,k] = deal(col{1},j);
    end    
    
    % turns the axis box on  
    hold(hAx{k},'on');
    set(hAx{k},'linewidth',1.5,'box','on')

    % sets the plot data values
    if isempty(p(j).I_mn)
        Yplt = NaN(size(Tplt));
        
    elseif isInact
        % case is inactivity proportion
        Yplt = p(j).I_mn;
        
    else
        % case is the hourly sleep duration
        Yplt = (1 - p(j).I_mn);
    end                
    
    % smooths the signal (if required)
    if pP.useSm
        Yplt = smooth(Yplt,pP.nSm);
    end

    % plots the signal SEM (if checked)
    if pP.plotErr                              
        % sets the error signal
        Yerr = p(j).I_sem(:);

        % sets the error signal and plots the error bars
        if ~isempty(p(j).I_mn)
            % plots the SEM signals
            tic
            set(gcf,'CurrentAxes',hAx{k})            
            plotSignalSEM(pY*Yplt,pY*Yerr,Tplt,colNw,0.5)
            toc
        end
    end

    % turns the grid on/off depending on the flag values
    if (j == 1) || ~sP.Sub.isComb
        if pP.plotGrid
            % grid is turned on
            grid on
        else
            % grid is turned off
            grid off
        end
        
        % sets the x/y axis limits
        if range(sP.xLim) == 0
            set(hAx{k},'ylim',[0 yLim],'box','off')            
        else
            set(hAx{k},'xlim',sP.xLim*tMlt,'ylim',[0 yLim],'box','off')
        end        
        
        % case is using the absolute time axis
        if pP.isZeitG         
            % case is using Zeitgeiber time
            setZeitGTimeAxis(hAx{k},p(j).T,snTot(1));
            
        else
            % if plotting day/night, then set absolute time
            setAbsTimeAxis(hAx{k},p(j).T,snTot(1));  
        end              
        
        % formats the plot axis
        formatPlotAxis(hAx{k},pF,i);
        set(hAx{k},'UserData',k);
        axis(hAx{k},'on')
    end
    
    % plots the traces (time scaled to hours)
    hPlotNw = plotFullSignal(hAx{k},Tplt,pY*Yplt);
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

% formats and resets the axis positions
resetAxesPos(hAx,m,n,[0 0.025]); 

% resets the legend location
if sP.Sub.isComb
    % creates the legend object
    hLg = createLegendObj(hPlot,pF.Legend);
    
    % resets the legend position
    pLg = get(hLg,'position');
    set(hLg,'Position',[(1-pLg(3)),0.5*(1-pLg(4)),pLg(3:4)])
    resetLegendPos(hLg,hAx)
end

% ----------------------------------------------------------------------- %
% ---                     MISCELLANEOUS FUNCTION                      --- %
% ----------------------------------------------------------------------- %

% --- sets up the time vectors over all experiments
function [snTot,Ttot] = setupTimeVectors(snTot)

% converts the time vectors
Ttot = arrayfun(@(x)(cell2mat(x.T)),snTot,'un',0);

% if there is more than one experiment, then reduce to the same length
if length(snTot) > 1
    Tmax = min(cellfun(@(x)(x(end)),Ttot));
    
    for i = 1:length(snTot)
        % determines the valid time points
        ii = Ttot{i} <= Tmax;
        
        % sets the feasible time/coordinate points
        Ttot{i} = Ttot{i}(ii);
        for j = find(~cellfun('isempty',snTot(i).iMov.flyok'))
            % resets the x-coordinates
            if ~isempty(snTot(i).Px)
                snTot(i).Px{j} = snTot(i).Px{j}(ii,:);
            end
            
            % resets the y-coordinates
            if ~isempty(snTot(i).Py)
                snTot(i).Py{j} = snTot(i).Py{j}(ii,:);
            end            
        end
    end
end

% --- interpolates the fly coordinates for the experiment
function snTot = interpFlyCoords(snTot,T,T0)

% interpolates the x-coordinates (if they exist)
if ~isempty(snTot.Px)
    % field retrieval
    for i = find(~cellfun('isempty',snTot.Px(:)'))
        snTot.Px{i} = cell2mat(cellfun(@(x)(...
            interp1(T0,x,T,'linear')),num2cell(snTot.Px{i},1),'un',0));
    end
end

% interpolates the x-coordinates (if they exist)
if ~isempty(snTot.Py)
    % field retrieval
    for i = find(~cellfun('isempty',snTot.Py(:)'))
        snTot.Py{i} = cell2mat(cellfun(@(x)(...
            interp1(T0,x,T,'linear')),num2cell(snTot.Py{i},1),'un',0));
    end
end

% --- sets up the time bin data arrays
function [T,indB,indD,useInt,ok] = setupTimeBinData(snTot,Ttot,cP,h)
            
% field retrieval and initialisations
ok = true;
useInt = false;
nExp = length(Ttot);
wOfs = nExp > 1;

% memory allocation
[T0,indB,indD] = deal(cell(nExp,1));

% sets up the time bin data for each experiment
for i = 1:nExp
    % updates the progressbar (if more than one expt)
    if wOfs
        wStr0 = 'Time Bin Setup';
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr0,i,nExp);        
        if h.Update(wOfs,wStrNw,i/(1+nExp))
            ok = false;
            return
        end        
    end
    
    % updates the progressbar
    h.Update(1+wOfs,'Calculating Time Bin Indices...',0);
    pause(0.01);    
    
    % determines the binned indices (for time length, tBin) and 
    % determines the bins which has at least two time points
    indB{i} = detTimeBinIndices(Ttot{i},cP.tMove);
    [T0{i},indB{i}] = setBinnedTimeArray(Ttot{i},indB{i},cP.tMove);
    
    % updates the progressbar
    h.Update(1+wOfs,'Setting Up Day/Time-Bin Mapping Indices...',0.5);
    pause(0.01);
    
    % determines daily time bin grouping indices (prop. inactivity only)
    iExpt = snTot(i).iExpt(1);
    indD{i} = detTimeGroupIndexArray(T0{i},iExpt,cP.tMove,cP.Tgrp0);
    
    % updates the progressbar
    h.Update(1+wOfs,'Time Bin Setup Complete!',1);
    pause(0.01);    
end

% sets the final time vectos
if nExp == 1
    % case is a single experiment
    T = T0{1};
    
else
    % case is multiple experiments
    
    % resets the day/time-bin mapping indices
    nD = min(cellfun(@(x)(size(x,2)),indD));
    indD = cellfun(@(x)(x(:,1:nD)),indD,'un',0);

    % resets the day/time-bin mapping indices to match over all expts
    indD0 = max(cellfun(@(x)(find(~isnan(x(:,1)),1,'first')),indD));
    indDF = min(cellfun(@(x)(find(~isnan(x(:,end)),1,'last')),indD));    
    for i = 1:nExp
        [indD{i}(1:(indD0-1),1),indD{i}((indDF+1):end,end)] = deal(NaN);
    end    
    
    % case is for multiple experiments
    if all(cellfun(@(x)(max(abs(T0{1}-x))),T0(2:end)) < cP.tMove)
        % case is all time vectors match
        T = T0{1};
        
    else
        % otherwise, reset the time vectors over all experiments
        iMn = argMin(cellfun('length',T0));
        [T,useInt] = deal(T0{iMn},true);        
        indD = repmat(indD(iMn),nExp,1);
    end
end

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)
