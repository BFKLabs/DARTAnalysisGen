% --- calculates the temporal proportional activity of the flies for each 
%     minute of the experiment
function pData = ProportionalActivity(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Temporal Activity Proportion';
pData.Type = {'Pop','Multi'};
pData.fType = [1 1 2 1];
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
    [pData.hasSP,pData.canComb] = deal(true);
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
nPara = 8;                      
pP = setParaFields(nPara);

% sets the parameter lists
pList = {'Activity','Inactivity'};
pList2 = {'Standard Error Mean','Standard Deviation'};
pList3 = {'Population','Individual'};

% parameter tab strings
a = {'1 - General','2 - Signal Smoothing'};

% sets the parameter fields
pP(1) = setParaFields(a{1},'List',{1,pList},'pMet','Plot Data Type');
pP(2) = setParaFields(a{2},'Boolean',0,'isSmooth','Plot Smoothed Signals');
pP(3) = setParaFields(a{2},'Number',5,'nSmooth','Signal Smooth Length',[1 15 true],{2,2});
pP(4) = setParaFields(a{1},'Number',1,'lWid','Plot Line Width',[0.1 10 false]);
pP(5) = setParaFields(a{1},'List',{1,pList3},'repType','Representation Type');
pP(6) = setParaFields(a{1},'List',{1,pList2},'errType','Error Type',[],{[5,7],[1,2]});
pP(7) = setParaFields(a{1},'Boolean',1,'plotErr','Plot Trace Errors',[],{5,1});
pP(8) = setParaFields(a{1},'Boolean',1,'plotGrid','Show Axis Gridlines');

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
oP = setupOutputParaStruct(snTot,false,true);
[Type1,Type2,xDep] = deal(5,4,{'Tbin'});

% sets the independent output variables
oP = addXVarField(oP,'Time Bin','Tbin','Time');
oP = addXVarField(oP,'Time','TT','Time');

% sets the dependent output variables
oP = addYVarField(oP,'Activity','Y',[],Type1,{'TT'},1);
oP = addYVarField(oP,'Activity (Mean)','Y_mn',[],Type2,xDep);
oP = addYVarField(oP,'Activity (SEM)','Y_sem',[],Type2,xDep);
oP = addYVarField(oP,'Activity (SD)','Y_sd',[],Type2,xDep);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
cP = retParaStruct(dcObj.pData.cP);
pP = retParaStruct(dcObj.pData.pP);
sP = retParaStruct(dcObj.pData.sP);

% sets the common class fields
dcObj.xName = 'Time';
dcObj.xUnits = 'min';
dcObj.yUnits = '%';
dcObj.yName = pP.pMet;
dcObj.useGrpHdr = true;
dcObj.combFig = sP.Sub.isComb;
dcObj.grpName = dcObj.pData.appName(sP.Sub.isPlot);

% sets up the class fields based on the representation type
if strcmp(pP.repType,'Population')
    % case is the population signals
    dcObj.pType = 'Trace';
    
else
    % case is the individual signals
    if dcObj.combFig
        dcObj.useGrpHdr = false;        
        dcObj.pType = 'Multi-Individual Trace';
    else
        dcObj.pType = 'Individual Trace';
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

% if using all the experiment for analysis, reset the calc parameters
if cP.useAll
    [cP.T0,cP.Tdur] = deal(0,floor(max(Tf)));
end
    
% sets the movement calculation type
cP.movType = 'Absolute Speed';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensions
nExp = length(snTot);

% memory allocation
nApp = length(snTot(1).iMov.ok);
plotD = initPlotValueStruct(snTot,pData,cP,...
                        'Tbin',[],'TT',cell(nExp,1),...
                        'Y_mn',[],'Y_sem',[],'Y_sd',[],'Y',[],'YT',[]);

% sets the group strings for each apparatus    
for i = 1:nApp
    plotD(i).Tbin = (0:(cP.Tdur-1))' + 0.5; 
    plotD(i).TT(:) = {plotD(i).Tbin};
end                               

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
    
    % calculates the video frame rate and experiment apparatus indices
    FPS = snTot(i).sgP.fRate/snTot(i).sgP.sRate;
    iApp = find(~cellfun('isempty',snTot(i).iMov.flyok));
    
    % sets the relevant time points and apparatus indices for this expt
    if cP.useAll
        % uses all the time points
        ii = 1:length(T{i});
    else
        % use only the points from the start to the duration end
        ii = (T{i} >= Tmlt*cP.T0) & (T{i} <= Tmlt*(cP.T0 + cP.Tdur));
    end    
    
    % sets the relevant time points and apparatus indices for this expt
    Tnw = T{i}(ii)-60*cP.T0;     
    isMove = calcFlyMove(snTot(i),Tnw,ii,iApp,cP.vAct);  
                    
    % calculates the time mid point of each frame, and from this determines
    % the time group that each frames belong to
    Tmid = (Tnw(1:end-1) + Tnw(2:end))/(2*Tmlt); 
    Tmid(diff(Tnw) > (2/FPS)*(60/Tmlt)) = -1;
    TmidF = floor(Tmid) + 1; 
    tGrp = cellfun(@(x)(find(TmidF == x)),num2cell(1:cP.Tdur)','un',0);                       
                
    % calculates the pre/post stimuli velocities for all flies, and
    % bins the values according to their time group bin
    for j = 1:length(iApp)  
        % calculates the speed values for the current group
        Z = cell2mat(cellfun(@(x)...
                            (sum(isMove{j}(x,:))/length(x)),tGrp,'un',0));

        % sets the values into the storage array
        if detMltTrkStatus(snTot(i).iMov)
            xiY = 1:size(Z,2);
            plotD(iApp(j)).Y(1,xiY,i) = num2cell(Z,1);
        else
            kk = 1:size(Z,2);
            plotD(iApp(j)).Y(1,kk,i) = num2cell(Z,1);            
            
%             fok = snTot(i).iMov.flyok{iApp(j)};
%             plotD(iApp(j)).Y(1,fok,i) = num2cell(Z,1);
        end
    end       
end

% sets/calculates the raw, mean and SEM proportional activity values
for i = 1:nApp
    % combines the group data
    plotD(i).YT = combGroupData(plotD(i).Y);
        
    % calculates the statistic metrics
    plotD(i) = calcMetricStats(plotD(i),'Y');    
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
cP = retParaStruct(pData.cP);
pP = retParaStruct(pData.pP);
sP = retParaStruct(pData.sP);
pF = pData.pF;

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the plotting indices and subplot indices
[ind,m,n,Tmlt] = deal(find(sP.Sub.isPlot),sP.Sub.nRow,sP.Sub.nCol,60);
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% sets the stimuli time stamps
T = cellfun(@(x)(cell2mat(x)),field2cell(snTot,'T'),'un',0);
Tf = cellfun(@(x)(x(end)),T)/Tmlt;

% if the experiment duration is small, then use second instead
if max(Tf) < 1
    [Tf,Tmlt] = deal(Tf*60,1); 
end

% converts the solution time arrays into single vectors
if cP.useAll
    [cP.T0,cP.Tdur] = deal(0,floor(max(Tf)));
end

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the subplot titles
if sP.Sub.isComb% || (all([m n] == 1))
    % case is combining, so remove the titles
    pF.Legend.String = snTot(1).iMov.pInfo.gName;
    for i = 1:length(pF.Title); pF.Title(i).String = ''; end
    
elseif all(cellfun('isempty',field2cell(pF.Title,'String')))
    % titles are empty, so reset them to the apparatus names
    for i = 1:length(pF.Title) 
        pF.Title(i).String = snTot(1).iMov.pInfo.gName{i}; 
    end
end

% retrieves the formatting struct
if isempty(m); szMx = 1; else; szMx = max([m n]); end
pF = retFormatStruct(pF,szMx);

% sets the legend strings (if combining onto a single figure
if sP.Sub.isComb
    [pF.Legend.String,m,n] = deal(pF.Legend.String(ind),1,1);
elseif ~all([m n] == 1)
    [pF.xLabel.ind,pF.yLabel.ind] = deal(NaN);  
end

% sets the x/y labels
[pF.xLabel.String,pF.yLabel.String] = deal('Time','% Activity');
pF.Title = pF.Title(ind);

% resets the x-label string
if Tmlt == 1
    pF.xLabel.String = sprintf('%s (sec)',pF.xLabel.String);
else
    pF.xLabel.String = sprintf('%s (min)',pF.xLabel.String);
end

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% memory allocation
hPlot = cell(nApp,1);
hAx = cell(1+(nApp-1)*(~sP.Sub.isComb),1);

% sets the trace colours
col = num2cell(distinguishable_colors(nApp,'w'),2);

% loops through all the subplots 
for j = 1:nApp
    % sets the subplot to be plotted on and updates the properties    
    if sP.Sub.isComb || (all([m n] == 1))
        % case is combining, so use the main axis
        colNw = col{j};
        if j == 1
            [hAx{j},hAxNw] = deal(createSubPlotAxes(hP)); hold on  
            axis(hAxNw,'on')
        end
    else
        % otherwise, create a seperate subplot          
        colNw = col{1};
        [hAx{j},hAxNw] = deal(createSubPlotAxes(hP,[m,n],j)); hold on  
    end    

    % sets the plot signals    
    if strcmp(pP.repType,'Population')
        if pP.isSmooth
            % data is smoothed
            Yplt = smooth(p(j).Y_mn*100,pP.nSmooth);
        else
            % data is not smoothed
            Yplt = p(j).Y_mn*100;
        end

        if strcmp(pP.pMet,'Inactivity'); Yplt = 100 - Yplt; end
        set(hAxNw,'linewidth',1.5,'box','on','UserData',j)

        % plots the SEM filled regions (if set)
        if pP.plotErr
            if strcmp(pP.errType,'Standard Error Mean')
                plotSignalSEM(Yplt,p(j).Y_sem*100,p(j).Tbin,colNw)        
            else
                plotSignalSEM(Yplt,p(j).Y_sd*100,p(j).Tbin,colNw)        
            end
        end

        % plots the traces     
        hPlot{j} = plot(hAxNw,p(j).Tbin,Yplt,'color',colNw,...
            'linewidth',pP.lWid,'UserData',j);
    else
        %
        Y0 = cellfun(@(x)(x.*(~isinf(x))),p(j).Y,'un',0); 
        Y1 = {reshape(Y0,[1 1 numel(Y0)])}; 
        Y2 = cellfun(@(x)(combineNumericCells3(x(~cellfun('isempty',x)))),Y1,'un',0); 
        
        % creates the individual plots
        Yplt = 100*squeeze(Y2{1});
        hPlot{j} = plot(hAxNw,p(j).Tbin,Yplt,'linewidth',pP.lWid);
        iPlt = arr2vec(1:length(hPlot{j}));
        
        if sP.Sub.isComb            
            arrayfun(@(h,i)(set...
                (h,'UserData',[j,i],'Color',col{j})),hPlot{j},iPlt);
        else
            arrayfun(@(h,i)(set(h,'UserData',i)),hPlot{j},iPlt);            
        end
    end

    % sets the x/y axis limits
    set(hAxNw,'xlim',[0 cP.Tdur]-0.01,'ylim',[0 100])        
    formatPlotAxis(hAxNw,pF,j);
    
    % adds in the gridlines (if checked)
    if pP.plotGrid; grid(hAxNw,'on'); end
end

% --- PLOT AXES REFORMATTING --- %
% ------------------------------ %

% sets the non-aligned x/y labels
if ~sP.Sub.isComb && (~all([m n] == 1))
    formatMultiXYLabels(hAx,pF,[m,n]);
else
    axis(hAxNw,'on')
end

% updates the figure if plotting a combined figure
if sP.Sub.isComb
    % resets the axis positions
    resetAxesPos(hAx,m,n,[0.02 0.02]);         
    
    % creates the legend object
    if sP.Sub.isComb && strcmp(pP.repType,'Individual')
        hPlotLg = cellfun(@(x)(x(1)),hPlot,'un',0);
        hLg = createLegendObj(hPlotLg,pF.Legend);
    else
        hLg = createLegendObj(hPlot,pF.Legend);
    end
    
    % resets the legend position           
    pLg = get(hLg,'position');
    set(hLg,'Position',[(1-pLg(3)),0.5*(1-pLg(4)),pLg(3:4)])   
    resetLegendPos(hLg,hAx)            
else
    % resets the axis positions
    resetAxesPos(hAx,m,n);     
end

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% converts the time bins from mins to secs
for i = 1:length(plotD)
    plotD(i).Tbin = plotD(i).Tbin*60;
end
