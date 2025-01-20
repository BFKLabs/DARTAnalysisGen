% --- plots the daily fly proportional activity/sleep mins per hour 
%    (long experiment only)
function pData = InactivityPlotDN(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Temporal Inactivity (Long)';
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
Tmin = 3;                       % short experiment max duration (in hours)
nPara = 4;
cP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the parameter list for 
pList = {'Hourly Sleep Duration','Inactivity Proportion'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'List',...
    {1,pList},'inactType','Inactivity Calculation Type');
cP(2) = setParaFields(a{1},'Number',...
    60,'tBin','Time Bin Size (min)',[1 120 false],{1,1});
cP(3) = setParaFields(a{1},'Number',...
    10,'tMove','Inactivity Detection Duration (sec)',[1 300 false],{1,2});

% sets the tool-tip strings
cP(1).TTstr = 'The method by which inactivity is calculated';
cP(2).TTstr = 'Time period over which sleep duration is calculated';
cP(3).TTstr = 'Time period over which inactivity proportion is calculated';

% if short experiments only, then create alignment checkbox option
T = field2cell(snTot,'T');
if all(convertTime(cellfun(@(x)(x{end}(end)),T),'sec','hrs') <= Tmin)
    % creates option
    cP(4) = setParaFields([],'Boolean',...
        0,'isAlign','Align All Experiments Start Points');
    cP(4).TTstr = ['Aligns all experiments to their respective ',...
                   'start times (Short experiments only)'];    
else
    % long experiment, so remove option
    cP = cP(1:3);
end

% --- initialises the plotting parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 5;
pP = setParaFields(nPara);

% sets the parameter lists
pList = {'SEM','Min/Max'};

% sets the tab list names
a = {'1 - General','2 - Signal Error'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},'Number',1,'lWid','Plot Line Width',[0.1 10 false]);
pP(2) = setParaFields(a{1},'Boolean',0,'isZeitG','Use Zeitgeiber Time Format');
pP(3) = setParaFields(a{1},'Boolean',0,'plotGrid','Show Plot Axis Gridlines');
pP(4) = setParaFields(a{2},'Boolean',0,'plotErr','Plot Signal Error');
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
dcObj.yName = cP.inactType;
dcObj.combFig = sP.Sub.isComb;
dcObj.tDay0 = [0,0,0,tDay,0,0];
dcObj.grpName = dcObj.pData.appName(sP.Sub.isPlot);

% sets the metric specific fields
switch cP.inactType
    case 'Hourly Sleep Duration'
        % case is the hourly sleep duration
        dcObj.yUnits = 'min/hour';
        
    case 'Inactivity Proportion'
        % case is the inactivity proportion
        dcObj.yUnits = '%';        
end

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
[nApp,nExp,ok] = deal(length(snTot(1).iMov.flyok),length(snTot),true);

% retrieves the parameter struct
[tBin,tMove] = deal(cP.tBin*60,cP.tMove);

% allocates memory for the temporary data/time plot arrays
if strcmp(cP.inactType,'Inactivity Proportion')
    tBinB = cP.tMove;
    T = ((tMove/2):tMove:convertTime(1,'day','sec'))';
else    
    tBinB = 60;
    T = ((tBin/2):(tBin):convertTime(1,'day','sec'))';
end

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,'T',T,'I',[],'I_mn',[],...
                             'I_sem',[],'I_min',[],'I_max',[],...
                             'inactType',cP.inactType,'hDay',gPara.TdayC);                           
                         
% creates the waitbar figure
wStr = {'Overall Progress','Inactivity Calculations'};
wOfs = (nExp > 1); wStr = wStr((2-wOfs):end);
h = ProgBar(wStr,'Inactivity Proportion Calculations');

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
        end
    end

    % determines the index offset
    Ttot = cell2mat(snTot(i).T);
    flyok = snTot(i).iMov.flyok;

    % determines the binned indices (for time length, tBin) and determines the
    % bins which has at least two time points
    h.Update(1+wOfs,'Determining Time Bin Indices...',0.25);
    indB = detTimeBinIndices(Ttot,tBinB);
            
    % calculates the binned activity 
    for j = 1:nApp
        % if there are no coordinates, then continue
        if isempty(snTot(i).Px{j})
            continue
        end
            
        % calculates the inactivity metrics based on the calculation type        
        if strcmp(cP.inactType,'Inactivity Proportion')
            % updates the waitbar figure
            wStrNw = sprintf(['Calculating Proportional Inactivity ',...
                              '(Region %i of %i)'],j,nApp);
            if h.Update(1+wOfs,wStrNw,0.50*(1+(j/(nApp+1))))
                % if the user cancelled, then exit the function
                [plotD,ok] = deal([],false);
                return                       
            end

            % determines the time group index array
            [Tnw,indB] = setBinnedTimeArray(Ttot,indB,cP.tMove);
            indD = detTimeGroupIndexArray(...
                            Tnw,snTot(i).iExpt(1),cP.tMove,cP.Tgrp0);

            % calculates the binned activity
            Anw = calcBinnedActivity(snTot(i),Ttot,indB,cP,j,flyok{j});    
            Inw = cellfun(@(x)(~x),Anw,'un',0);

            % sets the raw values into the storage array
            plotD(j) = setRawDataValues(...
                            plotD(j),snTot(i),Inw,indD,'I',i,j,1);

        else
            % updates the waitbar figure
            wStrNw = sprintf(['Calculating Sleep Metrics ',...
                              '(Region %i of %i)'],j,nApp);
            if h.Update(1+wOfs,wStrNw,0.50*(1+(j/(nApp+1))))
                % if the user cancelled, then exit the function
                [plotD,ok] = deal([],false);
                return                       
            end                    

            % sets the metric values
            cP2 = cP; 
            cP2.tMove = 5;
            [~,Inw] = calcSleepMetrics(...
                                snTot(i),Ttot,indB,cP2,j,flyok{j}); 

            % sets the raw values into the storage array
            plotD(j) = setRawDataValues(...
                                plotD(j),snTot(i),Inw,[],'I',i,j,2);
        end
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
    plotD(j) = calcMetricStats(plotD(j),'I');
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

% sets the time multiplier (from seconds to hours)
sP.xLim = [0 convertTime(24,'hrs','sec')];
[tMlt,absTime] = deal(convertTime(1,'sec','hrs'),true);

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

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

% sets the time axis properties
if strcmp(p(1).inactType,'Inactivity Proportion')
    [yLim,pY] = deal(100);
else
    % case is using Zeitgeiber time
    [yLim,pY] = deal(60,1);
    pF.yLabel.String = sprintf('Sleep (Min/Hour)'); 
end

% % updates the output parameter name and format struct
% pData.oP{4,1} = sprintf('%s (Mean)',p(1).inactType);
% pData.oP{5,1} = sprintf('%s (SEM)',p(1).inactType);
% pData.pF = pF;

% retrieves the formatting struct
if isempty(m); szMx = 1; else; szMx = max([m n]); end
pF = retFormatStruct(pF,szMx);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% memory allocation
hPlot = cell(nApp,1);

% other initialisations
Tax = [0 24];
dT = diff(p(1).T(1:2));
hAx = cell(1+(nApp-1)*(~sP.Sub.isComb),1);

% sets the trace colours
col = num2cell(distinguishable_colors(nApp,'w'),2);
    
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
        if (j == 1); plotSingleDNGraph(yLim,p(j)); end         
    else
        % otherwise, create a seperate subplot          
        colNw = col{1};
        [hAx{j},hAxNw] = deal(createSubPlotAxes(hP,[m,n],j)); hold on                        
        
        % adds the day/night (if not aligned)
        plotSingleDNGraph(yLim,p(j)); 
    end
    
    % turns the axis box on
    set(hAxNw,'linewidth',1.5,'box','on','UserData',j)

    % if there is an overlap, then plot the SEM traces
    Tplt = [(p(j).T(1)-dT/2);p(j).T;(p(j).T(end)+dT/2)]*tMlt;

    % sets the plot data values
    if isempty(p(j).I_mn)
        Yplt = NaN(size(Tplt));
    else
        YpltE = mean(p(j).I_mn([1 end]));
        Yplt = [YpltE;p(j).I_mn(:);YpltE];
    end

    % plots the signal SEM (if checked)
    if pP.plotErr                              
        switch pP.errType
            case 'SEM'
                % sets the error signal
                YpltS = reshape(p(j).I_sem,length(p(j).I_sem),1);            
            case 'Min/Max'
                % sets the error signal
                p(j).I_min = reshape(p(j).I_min,length(p(j).I_min),1);
                p(j).I_max = reshape(p(j).I_max,length(p(j).I_max),1);                                            
                YpltS = [p(j).I_min p(j).I_max];
        end

        % sets the error signal and plots the error bars
        if ~isempty(p(j).I_mn)
            Yend = mean(YpltS([1 end],:),1);               
            Yerr = [Yend;YpltS;Yend];       

            plotSignalSEM(pY*Yplt,pY*Yerr,Tplt,colNw,0.5)
        end
    end

    % plots the traces (time scaled to hours)   
    hPlotNw = plot(hAxNw,Tplt,pY*Yplt,'color',colNw,'linewidth',pP.lWid);
    if sP.Sub.isComb; hPlot{j} = hPlotNw; end
    set(hAxNw,'linewidth',0.5)        

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
        set(hAxNw,'xlim',[0 24],'ylim',[0 yLim])        
        
        % case is using the absolute time axis
        if pP.isZeitG         
            % case is using Zeitgeiber time
            setZeitGTimeAxis(hAxNw,convertTime(Tax,'hrs','sec'));        
        else
            % if plotting day/night, then set absolute time
            if absTime
                setAbsTimeAxis(hAxNw,convertTime(Tax,'hrs','sec'));  
            end
        end      

        % formats the plot axis
        formatPlotAxis(hAxNw,pF,i);
        expandAxesLim(hAxNw);
        axis(hAxNw,'on')
    end
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

% updates the figure if plotting a combined figure
if sP.Sub.isComb 
    % creates the legend object
    hLg = createLegendObj(hPlot,pF.Legend);
    
    % resets the legend position
    pLg = get(hLg,'position');
    set(hLg,'Position',[(1-pLg(3)),0.5*(1-pLg(4)),pLg(3:4)])
    resetLegendPos(hLg,hAx)
end

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)
    
% initialisations
cP = retParaStruct(pData.cP);

% re-orders the hourly sleep duration
if strcmp(cP.inactType,'Hourly Sleep Duration')
    % retrieves the variable string
    Var = field2cell(pData.oP.yVar,'Var');
    for i = 1:length(Var)
        if strcmp(Var{i},'I')
            for j = 1:length(plotD)
                Ynw = plotD(j).(Var{i});
                plotD(j).(Var{i}) = cellfun(@(x)(x(:)),Ynw,'un',0);
            end
        else
            for j = 1:length(plotD)
                plotD(j).(Var{i}) = arr2vec(plotD(j).(Var{i}));
            end
        end
    end
end