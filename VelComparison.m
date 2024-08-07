% --- comparison of the average pre/post-stimuli speeds for different time 
%     throughout the day (long experiment only)
function pData = VelComparison(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Pre & Post Stimuli Speed Comparison';
pData.Type = {'Pop','Multi'};
pData.fType = [2 2 3 1];
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
    [pData.hasSR,pData.hasRS] = deal(true,false);
    pData.sP = initSpecialPara(snTotL,pData,'nGrp',0);  
end
    
% ----------------------------------------------------------------------- %
% ---                 PARAMETER STRUCT SETUP FUNCTIONS                --- %
% ----------------------------------------------------------------------- %

% --- sets the function required information struct
function rI = initFuncReqInfo(pData)

% memory allocation
rI = struct('Scope',[],'Dur',[],'Shape',[],...
            'Stim',[],'Spec',[],'SpecFcn',[],'ClassicFcn',false);
        
% sets the struct fields
rI.Scope = setFuncScopeString(pData.Type);
rI.Dur = 'Long';
rI.Shape = 'None';
rI.Stim = 'Motor';
rI.Spec = 'None';        

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 3;                       
cP = setParaFields(nPara);
isLong = ~detIfShortExpt(field2cell(snTot,'T'));

% sets the tab list names
a = {'1 - General'};

% sets the parameter list for the drop-down list
pList = cellfun(@num2str,num2cell([1 2 4 6 8 12]),'un',0);

% sets the parameter fields
cP(1) = setParaFields(a{1},'List',{1+2*isLong,pList},'nGrp','Number of Daily Time Groups',[],{0,~isLong});
cP(2) = setParaFields(a{1},'Number',60,'tBefore','Pre-Stimuli Duration (sec)',[10 300 false]);
cP(3) = setParaFields(a{1},'Number',60,'tAfter','Post-Stimuli Duration (sec)',[10 300 false]);
cP(4) = setParaFields(a{1},'Boolean',0,'useZG','Use Zeitgeiber For Group Times');

% sets the other tool-tip strings
cP(1).TTstr = 'The number of groups that the day is split up into';
cP(2).TTstr = 'Time period before stimuli used to calculate average speed';
cP(3).TTstr = 'Time period after stimuli used to calculate average speed';

% adds the unique motor parameters
cP = addUniqueMotorPara(cP,snTot);

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 9;                       
pP = setParaFields(nPara);

% sets the parameter list strings
pList = {'Bar Graph','Boxplot'};

% sets the tab list names
a = {'1 - Scatterplot','2 - Metrics'};

% sets the parameter fields
pP(1) = setParaFields(a{2},'List',{1,pList},'pType','Plot Type');
pP(2) = setParaFields(a{1},'Number',6,'mSize','Plot Marker Size',[1 20 1]);
pP(3) = setParaFields(a{2},'Number',0.8,'pW','Bar Graph Relative Width',[0 1 false],{1,1});
pP(4) = setParaFields(a{2},'Boolean',0,'grpType','Group Data By Time Epoch');
pP(5) = setParaFields(a{2},'Boolean',1,'plotErr','Show Error Bars/Outliers');
pP(6) = setParaFields(a{1},'Boolean',0,'plotPolar','Plot Polar Coordinates');
pP(7) = setParaFields(a{1},'Boolean',1,'addTrend','Add Regression Line & Equation',[],{6,1});
pP(8) = setParaFields(a{2},'Boolean',0,'plotGrid','Plot Bar Graph/Boxplot Gridlines');
pP(9) = setParaFields(a{2},'Boolean',1,'plotFixedY','Fixed Subplot Limits To Overall Maximum');

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
[nApp,nSub] = deal(length(snTot.iMov.ok),3);
pF = setFormatFields(nSub);

% initialises the font structs
pF.Title = setFormatFields(setupFontStruct('FontSize',20),'',nSub);
pF.xLabel = setFormatFields(setupFontStruct('FontSize',14),'',nSub);
pF.yLabel = setFormatFields(setupFontStruct('FontSize',14),'',nSub);
pF.Axis = setFormatFields(setupFontStruct('FontSize',12),[]);

% sets x/y labels strings for the 1st plot
pF.xLabel(1).String = 'Pre-Stimuli Speed (mm s^{-1})';
pF.yLabel(1).String = 'Post-Stimuli Speed (mm s^{-1})';

% sets x/y labels strings for the 2nd plot
pF.yLabel(2).String = 'Speed (mm s^{-1})';

% sets x/y labels strings for the 2nd plot
pF.yLabel(3).String = 'Speed (mm s^{-1})';

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Axis(1).String{i} = snTot.iMov.pInfo.gName{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);
[Stats,Type] = deal({'CompMulti','Tgrp'},[1 3]);

% sets the independent output variables
oP = addXVarField(oP,'Time Group','Tgrp','Group');

% sets the dependent output variables
oP = addYVarField(oP,'Pre-Stim Spd','V1',Stats,Type,{'Tgrp'},1);
oP = addYVarField(oP,'Post-Stim Spd','V2',Stats,Type,{'Tgrp'},1);
oP = addYVarField(oP,'Speed Ratio','M',Stats,Type,{'Tgrp'},1);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(~,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% field retrievals
iAx = min(3,dcObj.getSelectAxesIndex);
pP = retParaStruct(dcObj.pData.pP);
sP = retParaStruct(dcObj.pData.sP);

% sets the common class fields
dcObj.grpName = dcObj.pData.appName;

% retrieves the currently selected axes index
if iAx == 1
    % case is the pre/post-stimuli speed scatterplot                
    
    % case is the bar graph plot types
    dcObj.pType = 'Scatterplot'; 
    dcObj.xName = 'Pre-Stimuli Speed';
    dcObj.xName2 = 'Time Group';
    dcObj.yName = 'Post-Stimuli Speed';    
    dcObj.xGrp = dcObj.plotD{1}(1).Tgrp(sP.pT);    
    [dcObj.xUnits,dcObj.yUnits] = deal('mm/sec');    
    
else
    % case is the bar graph plot types    
    
    % cell array fields
    Tgrp = dcObj.plotD{1}(1).Tgrp(sP.pT);
    yUnits = {'unitless','mm/sec'};
    pMet = {'Gradient','Average Speed'};
    
    % sets up the data cursor string type (based on plot type)
    if strcmp(pP.pType,'Bar Graph')
        % case is plotting bar graphs
        pType = {'Bar Graph','Multi-Bar Graph'};
    else
        % case is plotting boxplots
        pType = {'Boxplot','Multi-Boxplot'};        
    end
    
    % sets the data-cursor fields for the subplot
    dcObj.pType = pType{iAx-1};    
    dcObj.yName = pMet{iAx-1};
    dcObj.yUnits = yUnits{iAx-1};            
    [dcObj.xUnits,dcObj.yGrp] = deal([]);
    
    % sets the graph dependent fields
    if iAx == 3
        % case is the pre/post-stimuli average speeds
        
        % initialisations
        sType = {'Pre-Stimuli','Post-Stimuli'};
        [tStr,sStr] = deal('Time Group','Time Epoch');        
        
        % sets the primary/secondary x-value names/values
        if pP.grpType
            % case is the data is grouped by time-groups
            [dcObj.xName,dcObj.xName2] = deal(sStr,tStr);
            [dcObj.xGrp,dcObj.xGrp2] = deal(sType,Tgrp);
        else
            % case is the data is grouped by pre/post-stimuli
            [dcObj.xName,dcObj.xName2] = deal(tStr,sStr);            
            [dcObj.xGrp,dcObj.xGrp2] = deal(Tgrp,sType);
        end
    else
        % sets the x-value name/values
        [dcObj.xName,dcObj.xGrp] = deal('Time Group',Tgrp);
    end
end

% sets up the data cursor string
dTxt = dcObj.setupCursorString();

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

% sets the movement calculation type
cP.movType = 'Absolute Distance';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
nGrp = str2double(cP.nGrp);
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);

% fixed parameters
[tBefore,tAfter] = deal(cP.tBefore,cP.tAfter);

% sets the daily time group strings
Tgrp = setTimeGroupStrings(nGrp,cP.Tgrp0,cP.useZG); 

% retrieves the other calculation parameters (if they exist)
[devType,chType] = deal([]);
if isfield(cP,'devType'); devType = cP.devType; end
if isfield(cP,'chType'); chType = cP.chType; end

% memory allocation
plotD = initPlotValueStruct(snTot,pData,cP,NaN(1,nGrp),'Tgrp',Tgrp,...
                                'V1',[],'V1_mn',[],'V1_sem',[],...
                                'V2',[],'V2_mn',[],'V2_sem',[],...
                                'M',[],'M_mn',[],'M_sem',[]);                                     
                             
% ---------------------------------------------- %
% --- PRE/POST STIMULI VELOCITY CALCULATIONS --- %
% ---------------------------------------------- %

% sets the waitbar offset (is > 0 for more than one 
[wStr,wOfs] = deal({'Overall Progress'},1);
if nargin == 5
    % resets the strings
    [h,wOfs] = deal(varargin{1},wOfs+1);    
    wStr = [h.wStr(1),wStr];   
    
    % resets the waitbar figure fields
    for i = 2:length(wStr); h.Update(i,wStr{i},0); end

    % collapses the waitbar figure
    hh0 = getHandleSnapshot(findall(h));
    h.collapseProgBar(1);   
    
else
    % creates the waitbar figure
    h = ProgBar(wStr,'Pre/Post Stimuli Velocity Calculations');    
end

% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure (if more than one solution file)
    wStrNw = sprintf('%s (Experiment %i of %i)',wStr{wOfs},i,nExp);
    if h.Update(wOfs,wStrNw,i/(1+nExp))
        [plotD,ok] = deal([],false);
        return
    end
        
    % sets the video/stimuli time stamps into a single vector
    flyok = snTot(i).iMov.flyok;
    Ttot = cell2mat(snTot(i).T);
    Ts = getMotorFiringTimes(snTot(i).stimP,devType,chType);
    if ~isempty(Ts)
        % removes any stimuli that are less than the pre-stimuli phase
        Ts = Ts(Ts > cP.tBefore);
        
        % determines the indices of the stimuli events within the total
        % experiment, and determines what time groups that the stimuli
        % events took place in         
        T0 = snTot(i).iExpt(1).Timing.T0;
        indGrp = detTimeGroupIndices(Ts,T0,nGrp,cP.Tgrp0,true);
        
        % calculates the pre-stimuli velocities
        indV1 = cellfun(@(x)(find(Ttot>(x-tBefore),1,'first'):...
                    find(Ttot<x,1,'last')),num2cell(Ts),'un',0);
        indV2 = cellfun(@(x)(find(Ttot>x,1,'first'):find(...
                    Ttot<(x+tAfter),1,'last')),num2cell(Ts),'un',0);                    
                
        % calculates the pre/post stimuli velocities for all flies, and
        % bins the values according to their time group bin
        for j = 1:nApp
            % calculates the pre-stimuli velocities   
            if ~isempty(snTot(i).Px{j})
                V1 = calcBinnedFlyMovement(snTot(i),Ttot,indV1,cP,j,flyok{j},1);
                V1 = cell2mat(cellfun(@(x)(x/tBefore),V1,'un',0));
                                
                % calculates the post-stimuli velocities                    
                V2 = calcBinnedFlyMovement(snTot(i),Ttot,indV2,cP,j,flyok{j},1);
                V2 = cell2mat(cellfun(@(x)(x/tAfter),V2,'un',0));
                                
                % sets the raw data values
                plotD(j) = setRawDataValues(plotD(j),snTot(i),V1,indGrp,'V1',i,j,3);
                plotD(j) = setRawDataValues(plotD(j),snTot(i),V2,indGrp,'V2',i,j,3);                                                      
            end
        end
    end    
end

% sets the final pre/post stimuli event velocities (for each apparatus)
for i = 1:nApp            
    if ~all(cellfun('isempty',plotD(i).V1(:)))
        % sets/calculates the plot values
        plotD(i).M = cellfun(@(x,y)((x./y)),plotD(i).V2,plotD(i).V1,'un',0);
        
        % calculates the metric statistics for each of the values
        plotD(i) = calcMetricStats(plotD(i),'V1',3,false);        
        plotD(i) = calcMetricStats(plotD(i),'V2',3,false);
        plotD(i) = calcMetricStats(plotD(i),'M',3,false);                                  
    end    
end
    
% updates the waitbar figure
if ~h.Update(wOfs,'Pre/Post Stimuli Velocity Calculations Complete!',1)
    if nargin == 5
        resetHandleSnapshot(hh0)
    else
        h.closeProgBar()
    end
end

% ----------------------------------------------------------------------- %
% ---                        PLOTTING FUNCTION                        --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,varargout] = plotFunc(snTot,pData,plotD,ind)

% global variables
global regSz newSz

% retrieves the plotting paraeter struct
pP = retParaStruct(pData.pP);
sP = retParaStruct(pData.sP);
cP = retParaStruct(pData.cP);
pF = pData.pF;

% if the current data set is empty, then exit the loop
if all(isnan(plotD{1}(sP.pInd).V1_mn)); return; end

% sets any missing fields
if ~isfield(sP,'pT')
    nGrp = str2double(cP.nGrp);
    sP.pT = true(nGrp,1);
end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the indices to be plotted
iPlotT = find(sP.pT);
if isempty(iPlotT); return; end

% sets the parameters from the data structs
p = plotD{1}(sP.pInd);

% determines if the group count field and calculated values match
if length(p.Tgrp) ~= str2double(cP.nGrp)
    % if not then output an error message to screen
    mStr = sprintf(['Error! The calculated and selected daily ',...
                    'time groups do not match. You will need to ',...
                    're-run the function calculations.']); 
    
    
    % exits the function
    waitfor(msgbox(mStr,'Invalid Parameter Selection','modal'))
    return
end

% other parameters
hLg = [];
nGrp = length(iPlotT);
nGrpT = str2double(cP.nGrp);
[xi,Ymax] = deal(1:nGrp,500);
isBar = strcmp(pP.pType,'Bar Graph');
[c,mm,lWid,nTick,pWL] = deal('rk','+o*xsd^v><ph',2,6,0.85);

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the legend strings
[lStr1,lStr2] = deal(p.Tgrp,{'Pre-Stimuli','Post-Stimuli'});   
pF.Legend(1).String = lStr1(iPlotT);

% sets the main title string
pF.Title(1).String = pF.Axis(1).String{sP.pInd};
a = ~cellfun('isempty',p.V1) & ~cellfun('isempty',p.V2);

% sets the polar plot labels
if pP.plotPolar
    pP.addTrend = false;
    pF.yLabel(1).String = 'Distance (mm s^{-1})';
    pF.xLabel(1).String = 'Phase Angle (deg)';
end

% ------------------------------ %
% --- AVERAGE SPEED VS GROUP --- %
% ------------------------------ %

% sets up the sub-plot
hAx = createSubPlotAxes(hP,[2,2],4); 
hold(hAx,'on');
set(hAx,'Units','Normalized','box','on')

% creates the bar graph/boxplot figures
[hPlot,xTick] = plotDoubleBarBoxMetrics(hAx,p,{'V1','V2'},pP,iPlotT);

% creates the subplot legend
if pP.grpType
    % resets the legend/x-axis tick labels
    pF.Legend.String = lStr1(iPlotT);  
    if isBar; xTick = [1 2]; end           
    
    % updates the graph properties
    set(hAx,'xTick',xTick,'xticklabel',lStr2)               
else
    % resets the legend/x-axis tick labels
    [pF.Legend.String,pF.Legend.lgHorz] = deal(lStr2,true);      
    if isBar
        % case is the graph is a bar graph
        [xLim,xi] = deal([1 nGrp],1:nGrp);
    else
        % case is the graph is a boxplot graph
        [xLim,xi] = deal(xTick([1 end]) + 0.5*[-1 1],xTick);
        set(hAx,'xTick',xTick)
    end    
    
    % updates the graph properties
    set(hAx,'xlim',xLim + 0.5*[-1 1])        
end

% updates the axis properties
formatPlotAxis(hAx,pF,2); 

% creates the legend object
if ~pP.grpType || (nGrp > 1)
    hLg = createLegendObj(hPlot,pF.Legend);
end
    
% turns the grid on (if specified)
if pP.plotGrid; set(hAx,'ygrid','on'); end

% sets the 
if pP.plotFixedY
    % determines the maximum y-axis extent over all apparatus
    pp = plotD{1};    
    if isBar
        % retrieves and reduces down the data arrays
        V1mn = cellfun(@(x)(x(:,iPlotT)),field2cell(pp,'V1_mn'),'un',0);
        V1sem = cellfun(@(x)(x(:,iPlotT)),field2cell(pp,'V1_sem'),'un',0);
        V2mn = cellfun(@(x)(x(:,iPlotT)),field2cell(pp,'V2_mn'),'un',0);
        V2sem = cellfun(@(x)(x(:,iPlotT)),field2cell(pp,'V2_sem'),'un',0);        
        
        % determines the overall maxima
        V1mx = max(cellfun(@(x,y)(detOverallLimit(x+y)),V1mn,V1sem));
        V2mx = max(cellfun(@(x,y)(detOverallLimit(x+y)),V2mn,V2sem));
    else        
        % retrieves and reduces down the data arrays
        V1 = cellfun(@(x)(cellfun(@(y)...
            (y(:,iPlotT)),x,'un',0)),field2cell(pp,'V1'),'un',0);
        V2 = cellfun(@(x)(cellfun(@(y)...
            (y(:,iPlotT)),x,'un',0)),field2cell(pp,'V2'),'un',0);
        
        % determines the overall maxima
        V1mx = max(cellfun(@(x)(detOverallLimit(x)),V1));
        V2mx = max(cellfun(@(x)(detOverallLimit(x)),V2));
    end
                        
    % sets the axis limits for the two speed comparison subplots
    setStandardYAxis(hAx,[],nTick,detOverallLimit([V1mx V2mx]));
end

% resets the axes position
resetAxesPos(hAx,1,1);

% resets the location of the legend object
if ~isempty(hLg)
    [axP,lgP] = deal(get(hAx,'position'),get(hLg,'position'));
    if pP.grpType
        % case is that time is not being grouped
        lgP(2) = axP(2) + 0.5*(axP(4)-lgP(4));
        set(hLg,'position',lgP)
        
        %
        if isBar && pP.plotErr
            hErr = findall(hAx,'tag','hErr');
            cSzNw = hErr(1).CapSize*(axP(3)-pWL*lgP(3))/axP(3);
            
            for i = 1:length(hErr)
                hErr(i).CapSize = cSzNw;
            end
        end
        
        % resets the position of the axis
        resetObjPos(hAx,'width',-lgP(3),1)
    else
        % case is that time is being grouped
        lgP(1) = axP(1) + 0.5*(axP(3)-lgP(3));
        lgP(2) = sum(axP([2 4]));
        set(hLg,'position',lgP);        
    end
end

% adds the x-axis labels 
if ~pP.grpType
    setGroupString(hAx,pF,xi,lStr1(iPlotT),-90)
end

% --------------------------------------------- %
% --- POST/PRE-STIMULI SPEED RATIO VS GROUP --- %
% --------------------------------------------- %

% sets up the sub-plot
hAx2 = createSubPlotAxes(hP,[2,2],2); 
hold(hAx2,'on');

% creates the bar graph/boxplot
% yLim = plotBarBoxMetrics(hAx2,1:nGrp,p,'M',pP,[]);             
yLim = plotBarBoxMetrics(hAx2,iPlotT(:)',p,'M',pP,[]);     

% updates the axis properties
pF.yLabel(3).String = 'Gradient (Unitless)';
set(hAx2,'xtick',1:nGrp,'xticklabel',[],...
         'xLim',[1 nGrp]+0.5*[-1 1],'box','on')        
    
% formats the axis plot
formatPlotAxis(hAx2,pF,3); 

% turns the grid on (if specified)
if pP.plotGrid; set(hAx2,'ygrid','on'); end

% sets the 
if pP.plotFixedY
    % determines the maximum y-axis extent over all apparatus
    if isBar
        % retrieves and reduces down the data arrays
        Mmn = cellfun(@(x)(x(:,iPlotT)),field2cell(pp,'M_mn'),'un',0);
        Msem = cellfun(@(x)(x(:,iPlotT)),field2cell(pp,'M_sem'),'un',0);        
        
        % determines the overall maxima
        Mmx = max(cellfun(@(x,y)(detOverallLimit(x+y)),Mmn,Msem));
    else
        % retrieves and reduces down the data arrays
        M = cellfun(@(x)(cellfun(@(y)...
            (y(:,iPlotT)),x,'un',0)),field2cell(pp,'M'),'un',0);
        
        % determines the overall maxima                
        Mmx = max(cellfun(@(x)(detOverallLimit(x)),M));
    end                      
    
    % sets the axis limits for the two speed comparison subplots
    if Mmx > Ymax
        set(hAx2,'ylim',[0 Mmx],'yscale','log');
    else    
        setStandardYAxis(hAx2,[],nTick,Mmx);        
    end
else
    % determines if the y-axis scale is extremely large
    if max(get(hAx2,'ylim')) > Ymax
        % if so, then use a log scale
        set(hAx2,'yscale','log');
    else
        % otherwise, use a standard axis
        setStandardYAxis(hAx2,[],6,yLim);
    end    
end

% resets the axis dimensions
resetAxesPos(hAx2,1,1);
resetObjPos(hAx2,'bottom',0.015,1);
resetObjPos(hAx2,'height',-0.030,1);

% resets the axis dimensions so that match the last sub-plot
dL = retObjDimPos(hAx2,1,1) - retObjDimPos(hAx,1,1);
resetObjPos(hAx2,'left',retObjDimPos(hAx,1,1));
if pP.grpType
    resetObjPos(hAx2,'width',dL,1);
    setGroupString(hAx2,pF,xi,lStr1(iPlotT),-90);    
else
    resetObjPos(hAx2,'width',retObjDimPos(hAx,3,1));        
end

% pause for update...
pause(0.05);

% -------------------------------------- %
% --- PRE/POST STIMULI SPEED MARKERS --- %
% -------------------------------------- %

% memory allocations
hPlot = zeros(nGrp,1);

% creates a new subplot
hAxM = subplot(1,2,1,'Parent',hP,'UserData',1); 
hold(hAxM,'on');
axis(hAxM,'off');

%
V1 = cellfun(@(x)(mean(x,1,'omitnan')),p.V1,'un',0);
V2 = cellfun(@(x)(mean(x,1,'omitnan')),p.V2,'un',0);
[V1,V2] = deal(V1(~cellfun('isempty',V1(:))),V2(~cellfun('isempty',V2(:))));
[V1,V2] = deal(cell2mat(V1(:)),cell2mat(V2(:)));

% sets up the polar coordinate values (if plotting in polar coords)
if pP.plotPolar

end

% loops through all the time groups plotting the speeds
for i = 1:length(iPlotT)
    % determines the row/column index
    j = iPlotT(i);
    [ii,jj] = deal(mod(i-1,ceil(nGrp/2))+1,floor((i-1)/ceil(nGrp/2))+1);

    % plots the values (if any exist)
    if ~all(isnan(V1(:,j)))
        % sets up the x/y plot values
        if pP.plotPolar
            % calculates the distance/angles (first iteration only)
            if i == 1
                D = sqrt(V1.^2 + V2.^2);
                Phi = atan2(V1,V2)*(180/pi);
            end
            
            [xPlt,yPlt] = deal(Phi(:,j),D(:,j));
        else
            [xPlt,yPlt] = deal(V1(:,j),V2(:,j));
        end
        
        % creates the scatterplot
        cmStr = [c(jj),mm(ii)];
        hPlot(i) = plot(hAxM,xPlt,yPlt,cmStr,...
            'Markersize',pP.mSize,'UserData',i);
    end
end    

% updates the plot axis limits
if pP.plotPolar
    set(hAxM,'xlim',[0 90],'box','on')    
else
    % retrieves the pre/post-stimuli speeds
    [V1T,V2T] = field2cell(plotD{1},{'V1','V2'});    
    if pP.plotFixedY   
        % determines the max mean V1 values
        V1T = cellfun(@(x)(cellfun(@(y)...
                        (y(:,iPlotT)),x,'un',0)),V1T,'un',0);
        V1Tmn = cellfun(@(y)(cellfun(@(x)...
                        (mean(x,1,'omitnan')),y(:),'un',0)),V1T,'un',0);
        V1Tmn = cellfun(@(y)(max...
                        (max(cell2mat(y(~cellfun('isempty',y)))))),V1Tmn);
        
        % determines the max mean V2 values
        V2T = cellfun(@(x)(cellfun(@(y)...
                        (y(:,iPlotT)),x,'un',0)),V2T,'un',0);
        V2Tmn = cellfun(@(y)(cellfun(@(x)...
                        (mean(x,1,'omitnan')),y(:),'un',0)),V2T,'un',0);
        V2Tmn = cellfun(@(y)(max...
                        (max(cell2mat(y(~cellfun('isempty',y)))))),V2Tmn);
        
        % sets the overall limit 
        VLim0 = [detOverallLimit(V1Tmn),detOverallLimit(V2Tmn)];
        YmxNw = detOverallLimit(VLim0);
    else
        V1mx = max(arr2vec(V1(:,iPlotT)));
        V2mx = max(arr2vec(V2(:,iPlotT)));
        YmxNw = detOverallLimit([V1mx,V2mx]);
    end
    
    % sets the axis limits so that they are square
    set(hAxM,'xlim',[0 YmxNw],'ylim',[0 YmxNw],'box','on') 
end      

% turns the grid on (if specified)
if pP.plotGrid; grid(hAxM,'on'); end

% updates the axis properties
formatPlotAxis(hAxM,pF,1); 
axis(hAxM,'normal');

% resets the axis location
hPos = get(hAxM,'Position');
resetObjPos(hAxM,'left',-0.15*hPos(3),1)
resetObjPos(hAxM,'width',0.20*hPos(3),1)                
    
% plots the gradient lines (if specified)
if pP.addTrend
    % calculates the font-size ratio
    fR = max(0.8,min(newSz(3:4)./regSz(3:4))*get(0,'ScreenPixelsPerInch')/72);
    fSzT = setMinFontSize(pF.Axis(1).Font.FontSize*fR,'text');
    
    % calculates the linear fits to the lines
    [M_mn,CI,R2] = calcLinearFit(V1,V2);    
    xLim = get(hAxM,'xlim');

    % plots the values
    hPltT = [plot(hAxM,xLim,M_mn*xLim,'g','linewidth',lWid),...
             plot(hAxM,xLim,CI(1)*xLim,'g--','linewidth',lWid),...
             plot(hAxM,xLim,CI(2)*xLim,'g--','linewidth',lWid)];
    arrayfun(@(x)(set(x,'HitTest','off')),hPltT);  
    
    % resets the order of the marker/lines
    hPltM = setdiff(get(hAxM,'Children'),hPltT);
    set(hAxM,'Children',[hPltM(:);hPltT(:)])
    
    % sets the new text label
    nwText = [{sprintf('Gradient = %.2f %s %.2f',M_mn,'\pm',diff(CI)/2)};...
              {sprintf('(R^2 = %.2f)',R2)}];
    hText = text(0,0,nwText,'FontWeight',...
                 'bold','FontSize',fSzT,'HorizontalAlignment','center',...
                 'parent',hAxM,'tag','Other');
    
    % resets the position of the text object
    [tExt,dX] = deal(get(hText,'Extent'),0.025*diff(xLim));
    pText = [xLim(2)-(dX+tExt(3)/2),2*tExt(4)/3,0];
    set(hText,'position',pText,'Units','Normalized')    
end

% creates the legend
if nGrpT > 1
    isOut = sP.pT' & ~all(isnan(V1),1);
    [dY,lStrF] = deal(0.02,lStr1(isOut));
    hasPlt = hPlot > 0;

    if sum(isOut) == 1
        % only one valid time group
        hLg = legend(hPlot(hasPlt),lStrF,'location','NorthOutside');
        set(hLg,'box','off')
    else
        % more than one valid time group
        hLg = columnlegend(...
            hPlot(hasPlt),2,lStrF,'location','NorthOutside');
    end

    % resets the axis locations
    set(hLg,'tag','hLgScatter')
    [axPos,lgPos] = deal(get(hAxM,'position'),get(hLg,'position'));
        
    % resets the legend position
    lgPos(1) = axPos(1) + (axPos(3) - lgPos(3))/2;
    lgPos(2) = 1 - (dY + lgPos(4));
    set(hLg,'position',lgPos)
    
    % resets the axis height
    axPos(4) = (lgPos(2) + lgPos(4)/2) - (axPos(2) + 4*dY);
    set(hAxM,'position',axPos)
end

% sets the final axes properties
axis(hAxM,'on');
setObjVisibility(hAx,'on')

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)