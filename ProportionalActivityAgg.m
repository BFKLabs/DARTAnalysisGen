% --- calculates the mean proportional activity of the flies for each 
%     minute over the duration of the experiment (short experiment only)
function pData = ProportionalActivityAgg(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Aggregate Activity Proportion';
pData.Type = {'Pop','Multi'};
pData.fType = [1 1 2 1];
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
    [pData.hasSP,pData.hasRC] = deal(true,false);
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
rI.Dur = 'Short';
rI.Shape = 'None';
rI.Stim = 'None';
rI.Spec = 'None';        

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 4;                       
cP = setParaFields(nPara);

% sets the parameter fields
cP(1) = setParaFields([],'Number',2,'vAct','Activity Threshold (mm/s)',[0.1 10 false]);
cP(2) = setParaFields([],'Boolean',1,'useAll','Analyse Entire Experiment');
cP(3) = setParaFields([],'Number',0,'T0','Start Time (min)',[0 inf true],{2,1});
cP(4) = setParaFields([],'Number',30,'Tdur','Analysis Duration (min)',[1 inf true],{2,1});

% sets the tool-tip strings
cP(1).TTstr = 'The inter-frame velocity threshold used to indicate movement';
cP(2).TTstr = 'Determines whether or not to use the entire experiment for the analysis';
cP(3).TTstr = 'The time into the experiment where the analysis begins';
cP(4).TTstr = 'The duration period of the analysis';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 7;                      
pP = setParaFields(nPara);

% sets the parameter list strings

% sets the parameter lists
pList = {'Activity/Inactivity Percentage','Active/Inactive Minutes',...
          'Average/Active Speed','Displacement'};
pList2 = {'Both Metrics','First Metric','Second Metric'};
pList3 = {'Boxplot','Bar Graph'};

% parameter tab strings
a = {'1 - General','2 - Plot Metrics'};

% sets the parameter fields
pP(1) = setParaFields(a{1},'List',{1,pList2},'pDisp','Plot Metric');
pP(2) = setParaFields(a{1},'List',{1,pList3},'pType','Display Metric',[],{1,[1 2 3]});
pP(3) = setParaFields(a{2},'List',{1,pList},'pMet','Plot Type');
pP(4) = setParaFields(a{2},'Number',0.8,'pW','Bar Graph Relative Width',[0 1 false],{3,1});
pP(5) = setParaFields(a{2},'Boolean',1,'grpTime','Group Metrics By Type');
pP(6) = setParaFields(a{1},'Boolean',1,'plotGrid','Show Axis Gridlines');
pP(7) = setParaFields(a{1},'Boolean',1,'plotErr','Show Error Bars/Outliers');

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
pF = setFormatFields(1);

% initialises the font structs
pF.Title = setFormatFields([],'Aggregate Proportional Activity',1);
pF.xLabel = setFormatFields([],'');
pF.yLabel = setFormatFields([],'');
pF.Axis = setFormatFields([],[]);

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);
[Type1,Type2,Stats] = deal(3,2,{'Comp'});

% sets the dependent output variables
oP = addYVarField(oP,'% Active','A',Stats,Type1,[],1);
oP = addYVarField(oP,'% Active (Mean)','A_mn',[],Type2);
oP = addYVarField(oP,'% Active (SEM)','A_sem',[],Type2);
oP = addYVarField(oP,'% Active (Median)','A_md',[],Type2);
oP = addYVarField(oP,'% Active (Lower Quartile)','A_lq',[],Type2);
oP = addYVarField(oP,'% Active (Upper Quartile)','A_uq',[],Type2);
oP = addYVarField(oP,'Overall Avg. Speed','V',Stats,Type1,[],1);
oP = addYVarField(oP,'Overall Avg. Speed (Mean)','V_mn',[],Type2);
oP = addYVarField(oP,'Overall Avg. Speed (SEM)','V_sem',[],Type2);
oP = addYVarField(oP,'Overall Avg. Speed (Median)','V_md',[],Type2);
oP = addYVarField(oP,'Overall Avg. Speed (Lower Quartile)','V_lq',[],Type2);
oP = addYVarField(oP,'Overall Avg. Speed (Upper Quartile)','V_uq',[],Type2);
oP = addYVarField(oP,'Active Avg. Speed','Va',Stats,Type1,[],1);
oP = addYVarField(oP,'Active Avg. Speed (Mean)','Va_mn',[],Type2);
oP = addYVarField(oP,'Active Avg. Speed (SEM)','Va_sem',[],Type2);
oP = addYVarField(oP,'Active Avg. Speed (Median)','Va_md',[],Type2);
oP = addYVarField(oP,'Active Avg. Speed (Lower Quartile)','Va_lq',[],Type2);
oP = addYVarField(oP,'Active Avg. Speed (Upper Quartile)','Va_uq',[],Type2);
oP = addYVarField(oP,'% Inactive','I',Stats,Type1,[],1);
oP = addYVarField(oP,'% Inactive (Mean)','I_mn',[],Type2);
oP = addYVarField(oP,'% Inactive (SEM)','I_sem',[],Type2);
oP = addYVarField(oP,'% Inactive (Median)','I_md',[],Type2);
oP = addYVarField(oP,'% Inactive (Lower Quartile)','I_lq',[],Type2);
oP = addYVarField(oP,'% Inactive (Upper Quartile)','I_uq',[],Type2);
oP = addYVarField(oP,'Active Duration','Tw',Stats,Type1,[],1);
oP = addYVarField(oP,'Active Duration (Mean)','Tw_mn',[],Type2);
oP = addYVarField(oP,'Active Duration (SEM)','Tw_sem',[],Type2);
oP = addYVarField(oP,'Active Duration (Median)','Tw_md',[],Type2);
oP = addYVarField(oP,'Active Duration (Lower Quartile)','Tw_lq',[],Type2);
oP = addYVarField(oP,'Active Duration (Upper Quartile)','Tw_uq',[],Type2);
oP = addYVarField(oP,'Inactive Duration','Ts',Stats,Type1,[],1);
oP = addYVarField(oP,'Inactive Duration (Mean)','Ts_mn',[],Type2);
oP = addYVarField(oP,'Inactive Duration (SEM)','Ts_sem',[],Type2);
oP = addYVarField(oP,'Inactive Duration (Median)','Ts_md',[],Type2);
oP = addYVarField(oP,'Inactive Duration (Lower Quartile)','Ts_lq',[],Type2);
oP = addYVarField(oP,'Inactive Duration (Upper Quartile)','Ts_uq',[],Type2);
oP = addYVarField(oP,'Displacement','D',Stats,Type1,[],1);
oP = addYVarField(oP,'Displacement (Mean)','D_mn',[],Type2);
oP = addYVarField(oP,'Displacement (SEM)','D_sem',[],Type2);
oP = addYVarField(oP,'Displacement (Median)','D_md',[],Type2);
oP = addYVarField(oP,'Displacement (Lower Quartile)','D_lq',[],Type2);
oP = addYVarField(oP,'Displacement (Upper Quartile)','D_uq',[],Type2);

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

% converts the solution time arrays into single vectors
Tmlt = 60;
T = cellfun(@(x)(cell2mat(x)),field2cell(snTot,'T'),'un',0);
Tf = cellfun(@(x)(x(end)),T)/Tmlt;

% if the experiment duration is small, then use second instead
if max(Tf) < 1
    [Tf,Tmlt] = deal(Tf*60,1); 
end

% checks to see if the solution struct has the sub-region data struct
ok = checkFuncPara({'AnalysisTimeCheck'},cP,Tf);
if ~ok; plotD = []; return; end
    
% sets the movement calculation type
cP.movType = 'Absolute Speed';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensions
[nApp,nExp,ok] = deal(length(snTot(1).iMov.flyok),length(snTot),true);

% initialises the plot data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                            'A',[],'A_sem',[],'A_mn',[],...
                            'A_md',[],'A_lq',[],'A_uq',[],...
                            'V',[],'V_sem',[],'V_mn',[],...
                            'V_md',[],'V_lq',[],'V_uq',[],...
                            'Va',[],'Va_sem',[],'Va_mn',[],...
                            'Va_md',[],'Va_lq',[],'Va_uq',[],...                              
                            'I',[],'I_sem',[],'I_mn',[],...
                            'I_md',[],'I_lq',[],'I_uq',[],...
                            'Tw',[],'Tw_sem',[],'Tw_mn',[],...
                            'Tw_md',[],'Tw_lq',[],'Tw_uq',[],...
                            'Ts',[],'Ts_sem',[],'Ts_mn',[],...
                            'Ts_md',[],'Ts_lq',[],'Ts_uq',[],...
                            'D',[],'D_sem',[],'D_mn',[],...
                            'D_md',[],'D_lq',[],'D_uq',[]);                                             

% ------------------------------------------------------- %
% --- INTER-FRAME DISTANCE CALCULATION & THRESHOLDING --- %
% ------------------------------------------------------- %

% creates the waitbar figure
wStr = {'Overall Progress'};
h = ProgBar(wStr,'Proportial Activity Calculations');

% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure (if more than one solution file)
    if nExp > 1
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
        if h.Update(1,wStrNw,i/(1+nExp))
            [plotD,ok] = deal([],false);
            return
        end
    end
    
    % calculates the video frame rate
    FPS = snTot(i).sgP.fRate/snTot(i).sgP.sRate;
    iApp = find(~cellfun(@isempty,snTot(i).iMov.flyok));
        
    % sets the relevant time points and apparatus indices for this expt
    if cP.useAll
        % uses all the time points
        ii = 1:length(T{i});
    else
        % use only the points from the start to the duration end
        ii = (T{i} >= Tmlt*cP.T0) & (T{i} <= Tmlt*(cP.T0 + cP.Tdur));
    end
        
    % rearranges the time and x/y location arrays    
    Tnw = T{i}(ii)-Tmlt*cP.T0;
    [isMove,Vtot,Vact,Dact] = calcFlyMove(snTot(i),Tnw,ii,iApp,cP.vAct);
                            
    % calculates the time mid point of each frame, and from this determines
    % the (minute) time group that each frames belong to
    isOK = diff(Tnw) < (2/FPS)*(60/Tmlt);
    for j = 1:length(iApp)  
        % initialisations        
        [k,A] = deal(iApp(j),sum(isMove{j}(isOK,:),1)/sum(isOK));
        ifok = find(snTot(i).iMov.flyok{k});
        
        % sets the metric values
        plotD(k).Tw(1,ifok,i) = num2cell(Tmlt*A);
        plotD(k).Ts(1,ifok,i) = num2cell(Tmlt*(1-A));
        plotD(k).V(1,ifok,i) = num2cell(Vtot{j});
        plotD(k).Va(1,ifok,i) = num2cell(Vact{j});
        plotD(k).D(1,ifok,i) = num2cell(Dact{j}/1000);
        plotD(k).A(1,ifok,i) = num2cell(100*A);
        plotD(k).I(1,ifok,i) = num2cell(100*(1-A));
    end       
end

% calculates the metric statistics
for i = 1:nApp
    % calculates the statistic metrics
    plotD(i) = calcMetricStats(plotD(i),'A');
    plotD(i) = calcMetricStats(plotD(i),'I');
    plotD(i) = calcMetricStats(plotD(i),'Tw');
    plotD(i) = calcMetricStats(plotD(i),'Ts');
    plotD(i) = calcMetricStats(plotD(i),'V');
    plotD(i) = calcMetricStats(plotD(i),'Va');    
    plotD(i) = calcMetricStats(plotD(i),'D');  
end
    
% closes the waitbar figure
if ~h.Update(1,'Proportional Activity Calculations Complete!',1)
    if nargin < 5; h.closeProgBar(); end
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
ind = find(sP.Sub.isPlot);
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% parameters
yOfs = 0.01;

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% retrieves the formatting struct
pF = retFormatStruct(pF,1);

% sets the 
StrT = snTot(1).iMov.pInfo.gName(ind);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% memory allocation
hAx = createSubPlotAxes(hP); 
hold(hAx,'on')

% sets the plot values for the 
switch (pP.pMet)
    case ('Activity/Inactivity Percentage')       
        [yLim,StrM] = deal([0 100],{'Active','Inactive'});
        [pStr,pF.yLabel(1).String] = deal({'A','I'},'Percentage of Time');
        [tStrP,tStrS] = deal({'Activity','Inactivity'},'Percentage');        
    case ('Active/Inactive Minutes')
        [yLim,StrM] = deal([0 60],{'Active Time','Inactive Time'});
        [pStr,pF.yLabel(1).String] = deal({'Tw','Ts'},'Time (min/hour)');
        [tStrP,tStrS] = deal({'Activity','Inactivity'},'Duration');        
    case ('Average/Active Speed')
        [yLim,StrM] = deal(NaN,{'Average Speed','Active Speed'});
        [pStr,pF.yLabel(1).String] = deal({'V','Va'},'Speed (mm/sec)');
        [tStrP,tStrS] = deal({'Average','Active'},'Speed');       
    case ('Displacement')
        [yLim,StrM] = deal(NaN,{'Displacement'});
        [pStr,pF.yLabel(1).String] = deal({'D'},'Distance (m)');
        [tStrP,tStrS] = deal({[]},'Active Displacement');               
end

% sets the metric display indices
if (strcmp(pP.pMet,'Displacement'))
    [pInd,tStr] = deal(1,sprintf('%s %s',tStrP{1},tStrS));
else
    switch (pP.pDisp)
        case ('Both Metrics')
            pInd = [1 2];
            tStr = sprintf('%s/%s %s',tStrP{1},tStrP{2},tStrS);
        case ('First Metric')
            [pInd,tStr] = deal(1,sprintf('%s %s',tStrP{1},tStrS));
        case ('Second Metric')        
            [pInd,tStr] = deal(2,sprintf('%s %s',tStrP{2},tStrS));
    end
end
    
% sets the x-axis label and legend strings
if (pP.grpTime)
    % data is grouped by genotype
    [xStr,lStr] = deal(StrT,StrM');
else
    % data is grouped by metric
    [xStr,lStr] = deal(StrM(pInd)',StrT);
end

% creates the bar graph/boxplot figures
if length(pInd) == 1
    xTick = 1:nApp;
    plotBarBoxMetrics(hAx,xTick,p,pStr{pInd},pP,yLim);    
else
    [hPlot,xTick] = plotDoubleBarBoxMetrics(hAx,p,pStr(pInd),pP);
end
    
% removes all the text labels
delete(findall(hAx,'type','text'))

% sets the plot
set(hAx,'xtick',xTick);
if (isnan(yLim))
    yMax = max(get(hAx,'ylim'));
    set(hAx,'ylim',[0 yMax]); 
    setStandardYAxis(hAx(1),[],6,yMax);    
else
    set(hAx,'ylim',yLim); 
end
  
% adds in the gridlines (if checked)
if (pP.plotGrid); grid(hAx,'on'); end

% ------------------------------ %
% --- PLOT AXES REFORMATTING --- %
% ------------------------------ %

% sets the title string
pF.Title(1).String = tStr; 

% sets the x-axis labels
formatPlotAxis(hAx,pF,1);

% formats and resets the axis positions
resetAxesPos(hAx,1,1); 
if (pP.grpTime)
    setGroupString(hAx,pF,xTick,xStr,30);
else
    setGroupString(hAx,pF,xTick,xStr,0);
end

% creates the legend object
if (~pP.grpTime) || (length(pInd) > 1)
    % sets up and creates the legend
    pF.Legend.String = lStr;
    hLg = createLegendObj(hPlot,pF.Legend);
    
    % resets the legend/axis position
    [lgPos,axPos] = deal(get(hLg,'position'),get(hAx,'position'));
    set(hLg,'position',[1-(lgPos(3)+0.01),0.5-(lgPos(4)/2),lgPos(3:4)])
    set(hAx,'position',[axPos(1:2),axPos(3)-(lgPos(3)+2*yOfs),axPos(4)])    
end

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)
