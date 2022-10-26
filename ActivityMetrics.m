% --- calculates the fly activity initiation (number of starts per second 
%     stopped) over the duration of an experiment (short experiment only)
function pData = ActivityMetrics(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Activity Metrics';
pData.Type = {'Pop','Multi'};
pData.fType = [1 1 2 1];
pData.rI = initFuncReqInfo(pData);

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
nPara = 7;                       
cP = setParaFields(nPara);
is2D = all(cellfun(@(x)(is2DCheck(x.iMov)),num2cell(snTot)));

% sets the tab list names
a = {'1 - General','2 - Thresholds'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'Boolean',0,'useRegion','Split Activity By Region',[],{0,~is2D});
cP(2) = setParaFields(a{1},'Boolean',1,'useAll','Analyse Entire Experiment');
cP(3) = setParaFields(a{1},'Number',0,'T0','Start Time (min)',[0 inf true],{2,1});
cP(4) = setParaFields(a{1},'Number',60,'Tdur','Analysis Duration (min)',[1 inf true],{2,1});
cP(5) = setParaFields(a{2},'Number',3,'vAct','Activity Threshold (mm/s)',[0.1 10 false]);
cP(6) = setParaFields(a{2},'Number',0.2,'tMove','Movement Start Duration (s)',[0.1 10 false]);
cP(7) = setParaFields(a{2},'Number',3,'rTol','Outer Edge Region Distance (mm)',[1 10 false],{1,2});

% sets the tool-tip strings
cP(1).TTstr = 'Determines whether the activity is separated by region';
cP(2).TTstr = 'Determines whether or not to use the entire experiment for the analysis';
cP(3).TTstr = 'The time into the experiment where the analysis begins';
cP(4).TTstr = 'The duration period of the analysis';
cP(5).TTstr = 'The inter-frame velocity threshold used to indicate movement';
cP(6).TTstr = 'The duration the fly has to move to be considered a start';
cP(7).TTstr = 'Distance the outer edge is from the circle circumference';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 5;
pP = setParaFields(nPara);

% sets the parameter list strings
pList = {'Action Initiation','Inter-Bout Interval','Mean Bout Length'};
pList2 = {'Automatic','Linear','Log'};

% sets the parameter fields
pP(1) = setParaFields([],'List',{1,pList},'pMet','Activity Metric');
pP(2) = setParaFields([],'List',{1,pList2},'yScale','Y-Axis Scale');
pP(3) = setParaFields([],'Boolean',1,'plotGrid','Plot Trace Gridlines');
pP(4) = setParaFields([],'Boolean',1,'plotErr','Plot Boxplot Outliers');
pP(5) = setParaFields([],'Boolean',1,'grpType','Group Metrics By Type');

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
pF = setFormatFields(1);

% initialises the font structs
pF.Title = setFormatFields([],'',1);
pF.xLabel = setFormatFields([],'');
pF.yLabel = setFormatFields([],'');
pF.Axis = setFormatFields([],[]);

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot,false,true);
[xDep1,xDep2,Type1,Type2] = deal({'ioStr'},{'iotStr'},3,2);
[Stats1,Stats2] = deal({'CompMulti','ioStr'},{'CompMulti','iotStr'});

% sets the independent variable fields
oP = addXVarField(oP,'Activity Region','ioStr','Group');
oP = addXVarField(oP,'Activity Region','iotStr','Group');

% sets the dependent output variables
oP = addYVarField(oP,'Start Count','NS',Stats1,Type1,xDep1,1);
oP = addYVarField(oP,'Stop Time','TS',Stats1,Type1,xDep1,1);
oP = addYVarField(oP,'Activity Initiation','AI',Stats1,Type1,xDep1,1);
oP = addYVarField(oP,'Activity Initiation (Mean)','AI_mn',[],Type2,xDep1);
oP = addYVarField(oP,'Activity Initiation (SEM)','AI_sem',[],Type2,xDep1);
oP = addYVarField(oP,'Activity Initiation (Median)','AI_md',[],Type2,xDep1);
oP = addYVarField(oP,'Activity Initiation (Lower Quartile)','AI_lq',[],Type2,xDep1);
oP = addYVarField(oP,'Activity Initiation (Upper Quartile)','AI_uq',[],Type2,xDep1);
oP = addYVarField(oP,'Mean Bout Length','MBL',Stats2,Type1,xDep2,1);
oP = addYVarField(oP,'Mean Bout Length (Mean)','MBL_mn',[],Type2,xDep2);
oP = addYVarField(oP,'Mean Bout Length (SEM)','MBL_sem',[],Type2,xDep2);
oP = addYVarField(oP,'Mean Bout Length (Median)','MBL_md',[],Type2,xDep2);
oP = addYVarField(oP,'Mean Bout Length (Lower Quartile)','MBL_lq',[],Type2,xDep2);
oP = addYVarField(oP,'Mean Bout Length (Upper Quartile)','MBL_uq',[],Type2,xDep2);
oP = addYVarField(oP,'Inter-Bout Interval','IBI',Stats1,Type1,xDep1,1);
oP = addYVarField(oP,'Inter-Bout Interval (Mean)','IBI_mn',[],Type2,xDep1);
oP = addYVarField(oP,'Inter-Bout Interval (SEM)','IBI_sem',[],Type2,xDep1);
oP = addYVarField(oP,'Inter-Bout Interval (Median)','IBI_md',[],Type2,xDep1);
oP = addYVarField(oP,'Inter-Bout Interval (Lower Quartile)','IBI_lq',[],Type2,xDep1);
oP = addYVarField(oP,'Inter-Bout Interval (Upper Quartile)','IBI_uq',[],Type2,xDep1);
          
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
T = cellfun(@(x)(cell2mat(x)),field2cell(snTot,'T'),'un',0);
Tf = cellfun(@(x)(x(end)),T)/60;

% checks to see if the solution struct has the sub-region data struct
ok = checkFuncPara({'AnalysisTimeCheck'},cP,Tf);
if (~ok); plotD = []; return; end

% sets the movement calculation type
cP.movType = 'Absolute Speed';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
nR = cP.useRegion + 1;
[nApp,nExp,ok] = deal(length(snTot(1).iMov.flyok),length(snTot),true);

% sets the region strings
ioStr = {'Inner','Outer'};
iotStr = {'Inner','Outer','Transition'};

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                            'AI',[],'AI_mn',[],'AI_sem',[],...
                            'AI_md',[],'AI_lq',[],'AI_uq',[],...
                            'MBL',[],'MBL_mn',[],'MBL_sem',[],...
                            'MBL_md',[],'MBL_lq',[],'MBL_uq',[],...
                            'IBI',[],'IBI_mn',[],'IBI_sem',[],...
                            'IBI_md',[],'IBI_lq',[],'IBI_uq',[],...
                            'NS',[],'TS',[],'ioStr',ioStr,'iotStr',iotStr);

% other memory allocations                             
[I,OE] = deal(cell(nApp,nExp));
TT = cell(1,nExp);
indMBL = [1 3 2];

% creates the waitbar figure
wStr = {'Movement Detection','Activity Metrics Calculations'};
h = ProgBar(wStr,'Activity Calculations');

% ------------------------------------------------------- %
% --- INTER-FRAME DISTANCE CALCULATION & THRESHOLDING --- %
% ------------------------------------------------------- %
    
% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure 
    wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
    if h.Update(1,wStrNw,i/(1+nExp))
        [plotD,ok] = deal([],false);
        return
    end
    
    % calculates the video frame rate and experiment apparatus indices
    FPS = snTot(i).sgP.fRate/snTot(i).sgP.sRate;
    iApp = find(~cellfun(@isempty,snTot(i).iMov.flyok));
    
    % sets the relevant time points and apparatus indices for this expt
    if (cP.useAll)
        % uses all the time points
        ii = 1:length(T{i});
    else
        % use only the points from the start to the duration end
        ii = (T{i} >= 60*cP.T0) & (T{i} <= 60*(cP.T0 + cP.Tdur));
    end    
    
    % sets the relevant time points and apparatus indices for this expt
    Tnw = T{i}(ii)-60*cP.T0;     
    isMove = calcFlyMove(snTot(i),Tnw,ii,iApp,cP.vAct);      
        
    % determines points where flies are on the outer edge (if required)
    if cP.useRegion
        % retrieves the 2D coordinates wrt to the circle centre
        [dPx,dPy,R] = get2DCoordsBG(snTot(i),iApp);
        if (length(iApp) == 1); [dPx,dPy,R] = deal({dPx},{dPy},{R}); end
                    
        % determines the edge position flags for each fly (over all frames)
        [sFac,jj] = deal(snTot(i).sgP.sFac,ii(1:end-1));
        onEdge = cellfun(@(x,y,r)(detFlyEdgePos(...
                            x(jj,:),y(jj,:),r,cP,sFac)),dPx,dPy,R,'un',0);
    end
    
    % calculates the time mid point of each frame, and removes all the
    % frames where the inter-frame time difference is large (ie, the
    % change-over in video recordings)
    [Tmid,isOK] = deal((Tnw(1:end-1) + Tnw(2:end))/2,diff(Tnw) < 2/FPS);
    
    % sets the 
    TT{i} = Tmid(isOK);
    for j = 1:length(isMove)
        % sets the movement flags
        I{iApp(j),i} = isMove{j}(isOK,:);
        I{iApp(j),i}(1,:) = false;
        
        % removes any short duration movement events
        for k = 1:size(I{iApp(j),i},2)
            iGrp = getGroupIndex(I{iApp(j),i}(:,k));
            ii = cellfun(@(x)((Tnw(x(end))-Tnw(x(1)-1)) <= cP.tMove),iGrp);
            I{iApp(j),i}(cell2mat(iGrp(ii)),k) = false;
        end
        
        % sets the outer-edge flags (if required)
        if (cP.useRegion); OE{iApp(j),i} = onEdge{j}(isOK,:); end
    end
end

% updates the waitbar figure
h.Update(1,sprintf('%s (Complete)',wStr{1}),1);
    
% sets/calculates the raw, mean and SEM proportional activity values
for i = 1:nApp
    % updates the waitbar figure (if more than one solution file)
    wStrNw = sprintf('%s (Apparatus %i of %i)',wStr{2},i,nApp);
    if h.Update(2,wStrNw,i/(1+nApp))
        [plotD,ok] = deal([],false);
        return
    end    
            
    % sets all of the experiments into a single array
    for j = 1:nExp        
        if ~isempty(I{i,j})
            % sets the acceptance flags            
            ifok = find(snTot(j).iMov.flyok{i});
            
            % groups the frames where the fly has moved. from this,
            % calculate the duration of these movement events.
            moveGrp = cellfun(@(x)(getGroupIndex(x)),...
                            num2cell(I{i,j},1),'un',0);
            stopGrp = cellfun(@(x)(getGroupIndex(~x)),...
                            num2cell(I{i,j},1),'un',0);          
                        
            % removes the first group (remove stationary flies and the
            % first frame which is counted as being an inactive frame)
            stopGrp = cellfun(@(x)(x(2:end)),stopGrp,'un',0);
            for k = 1:length(moveGrp)
                if (~isempty(moveGrp{k}))
                    if (moveGrp{k}{end}(end) == size(I{i,j},1))
                        moveGrp{k} = moveGrp{k}(1:end-1);
                    end
                end
            end
            
            % calculates the duration of the stopped events
            tGrpStop = cellfun(@(x)(cellfun(@(y)(TT{j}(y(end)) - ...
                            TT{j}(y(1)-1)),x)),stopGrp,'un',0);            
            kk = find(~cellfun(@isempty,tGrpStop));     
            
            % calculates the duration of the moved events
            tGrpMove = cellfun(@(x)(cellfun(@(y)(TT{j}(y(end)) - ...
                            TT{j}(y(1)-1)),x)),moveGrp,'un',0);                                  
            jj = find(~cellfun(@isempty,tGrpMove));            
            
            % sets the region indices for each stop group
            [iIBI,iMBL,iAI] = deal(cell(size(stopGrp)));
            if (cP.useRegion)
                % determines the locations of the stopped events
                iIBI(kk) = cellfun(@(x,y)(cellfun(@(xx)(mean(x(xx)) > 0.5),y)+1),...
                            num2cell(OE{i,j}(:,kk),1),stopGrp(kk),'un',0);
                        
                % determines the location of the start events. from this, 
                % set the MBI/AI indices
                a = cellfun(@(x,y)(cell2mat(cellfun(@(xx)(...
                            x([xx(1),xx(end)]+[-1 1])'),y(:),'un',0))),...
                            num2cell(OE{i,j}(:,jj),1),moveGrp(jj),'un',0);     
                iMBL(jj) = cellfun(@(x)(sum(x,2)+1),a,'un',0);
                iAI(jj) = cellfun(@(x)(x(:,1)+1),a,'un',0);
            else
                % data is not separated by regions, so return ones array                 
                iIBI(kk) = cellfun(@(x)(ones(size(x))),stopGrp(kk),'un',0);
                [iMBL(jj),iAI(jj)] = deal(...
                          cellfun(@(x)(ones(size(x))),moveGrp(jj),'un',0));
            end                               
            
            % calculates the inter-bout intervals
            for k = 1:nR
                % determines the stopped events that correspond to the
                % region of interest
                ii = cellfun(@(x)(x == k),iIBI,'un',0);
                
                % sets the raw/mean IBI values for each fly in the current
                % experiment (only if fly is within specific region)
                IBI = cellfun(@(x,y)(mean(x(y),'omitnan')),tGrpStop,ii);
                for k2 = 1:length(ifok)
                    plotD(i).IBI{1,ifok(k2),j}(k) = IBI(k2);
                end
            end                                                                                     
                                    
            % calculates the inter-bout intervals            
            for k = 1:(1+2*(nR==2))
                % determines the start events that correspond to the
                % region of interest
                ii = cellfun(@(x)(x == k),iMBL,'un',0);            
            
                % sets the mean MBL values for each fly in the current
                % experiment (only if fly is within specific region)
                MBL = cellfun(@(x,y)(mean(x(y),'omitnan')),tGrpMove,ii);
                for k2 = 1:length(ifok)
                    plotD(i).MBL{1,ifok(k2),j}(indMBL(k)) = MBL(k2);
                end
            end
              
            % calculates the action initiation values for each fly
            for k = 1:nR         
                % determines the start events that correspond to the
                % region of interest, and determines their count
                ii = cellfun(@(x)(x == k),iAI,'un',0);       
                nStart = cellfun(@sum,ii);
            
                % sets the mean AI values for each fly in the current
                % experiment (only if fly is within specific region)
                [tStop,dnS] = deal(NaN(size(nStart)),zeros(size(nStart)));
                for j2 = jj
                    if (nStart(j2) > 0)
                        % removes the start event (if it is the first)
                        ii2 = find(ii{j2});
                        if (ii2(1) == 1)
                            dnS(j2) = 1;
                            [nStart(j2),ii2] = deal(nStart(j2)-1,ii2(2:end));
                        end
                        
                        % calculates the total stopped time for the events
                        % preceding the start events
                        if (~isempty(ii2))
                            tStop(j2) = sum(tGrpStop{j2}(ii2-1));
                        end
                    end
                end    
                
                % calculates the final AI values
                for k2 = 1:length(ifok)
                    plotD(i).AI{1,ifok(k2),j}(k) = nStart(k2)/tStop(k2);
                    plotD(i).NS{1,ifok(k2),j}(k) = nStart(k2)+dnS(k2);
                    plotD(i).TS{1,ifok(k2),j}(k) = tStop(k2);
                end
            end                    
        end        
    end
    
    % calculates the metric statistics
    plotD(i) = calcMetricStats(plotD(i),'AI');
    plotD(i) = calcMetricStats(plotD(i),'MBL');
    plotD(i) = calcMetricStats(plotD(i),'IBI');
end
    
% closes the waitbar figure
if ~h.Update(2,'Activity Calculations Complete!',1)
    if nargin < 5; h.closeProgBar(); end
end

% ----------------------------------------------------------------------- %
% ---                        PLOTTING FUNCTION                        --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,varargout] = plotFunc(snTot,pData,plotD,ind)

% retrieves the plotting parameter struct
pP = retParaStruct(pData.pP);
sP = retParaStruct(pData.sP);
cP = retParaStruct(pData.cP);
pF = pData.pF;

% sets the stimuli index
hPara = findall(0,'tag','figAnalysisPara');
hCheck = findall(hPara,'tag','grpType');
setObjEnable(hCheck,cP.useRegion);

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% turns off all the warnings
wState = warning('off','all');

% sets the plotting indices and subplot indices
ind = find(sP.Sub.isPlot);
nApp = length(ind); if (nApp == 0); return; end
[p,addLg] = deal(plotD{1}(ind),false);

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% retrieves the formatting struct
pF = retFormatStruct(pF,1);
xi = 1:length(ind);

% sets the title string
pF.Title.String = pP.pMet;

% sets the x-axis strings
if isfield(snTot(1).iMov,'Name')
    xStr = snTot(1).iMov.Name(ind);
else
    xStr = snTot(1).iMov.pInfo.gName(ind);
end

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% memory allocation
hAx = createSubPlotAxes();
% set(hAx,'Units',get(get(hAx,'parent'),'Units'))

% sets the plot values for the 
switch (pP.pMet)
    case ('Action Initiation')
        [lStr,pStr] = deal({'Inner','Outer'},'AI');
        pF.yLabel.String = 'Initiation (Starts/sec)';
    case ('Inter-Bout Interval')
        [lStr,pStr] = deal({'Inner','Outer'},'IBI');
        pF.yLabel.String = 'Duration (sec)';
    case ('Mean Bout Length')        
        [lStr,pStr] = deal({'Inner','Outer','Transition'},'MBL');
        pF.yLabel.String = 'Duration (sec)';
end

% creates the box-plot
if cP.useRegion
    % case is for specific region metrics
    pP.pType = 'Boxplot';
    [hPlot,xi,Yplt] = plotMultiBarBoxMetrics(hAx,p,pStr,pP); 
    yL0 = [min(Yplt(:)),max(Yplt(:))];
    
    % sets the legend strings
    if (pP.grpType)
        pF.Legend.String = lStr;
    else
        pF.Legend.String = snTot(1).iMov.pInfo.gName;
        xStr = lStr;
    end
    
    % creates the legend objects
    addLg = length(hPlot) > 1;
else
    % case is for specific non-region metrics
    Yplt = combineNumericCells(field2cell(p,pStr));  
    
    % creates the box plot
    if (pP.plotErr)
        % outliers are included
        boxplot(Yplt,'sym','r*');    
        yL0 = [min(Yplt(:)),max(Yplt(:))];
    else
        % outliers are not included
        hPlot = boxplot(Yplt,'sym','r');  
        yL0 = getBoxPlotLimits(hPlot,false);
        set(hAx,'ylim',[min(get(hAx,'ylim')),ceil(yL0(2))])
    end             
    
    % removes all the text labels
    set(hAx,'xlim',[0,size(Yplt,2)]+0.5)
    delete(findall(hAx,'type','text'))        
end
        
% sets the plot
if (strcmp(pP.yScale,'Automatic'))
    yScaleStr = {'linear','log'};
    yScale = yScaleStr{1 + (diff(log10(yL0))>2)};
else
    yScale = lower(pP.yScale);
end

%
set(hAx,'xtick',xi,'yscale',yScale,'UserData',1)
if (strcmp(yScale,'log'))
    yL = [10^floor(log10(yL0(1))) 10^ceil(log10(yL0(2)))];    
    pF.yLabel.String = sprintf('log_{10} %s',pF.yLabel.String);
else 
    yL = [0 max(get(hAx,'ylim'))];
end

% adds in the gridlines (if checked)
set(hAx,'ylim',yL);
if (pP.plotGrid) 
    set(hAx,'ygrid','on'); 
    xG = 0.5*(xi(1:end-1)+xi(2:end));

    for i = 1:length(xG)
        plot(xG(i)*[1 1],yL,'k:')
    end
end  
      
% ------------------------------ %
% --- PLOT AXES REFORMATTING --- %
% ------------------------------ %

% sets the x-axis labels
formatPlotAxis(hAx,pF,1);

% formats and resets the axis positions
resetAxesPos(hAx,1,1); 

% sets the group strings
setGroupString(hAx,pF,xi,xStr,30);

% resets the legend position (if one is created)
if (addLg)
    %
    hLg = createLegendObj(hPlot,pF.Legend,1,0);          
    
    % creates the legend object 
    [lgP,dX] = deal(resetVertLegendWidth(hLg),0.0);
    lgP(1:2) = [(1-lgP(3)+dX),(0.5-lgP(4)/2)];
    set(hLg,'position',lgP);      

    % resets the width of the trace plot axis
    pAx = get(hAx,'position');
    lgPNw = get(hLg,'position');
    resetObjPos(hAx,'Width',lgPNw(1)-(pAx(1)+dX))
end
    
% reverts the warning back to their original state
warning(wState);

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% retrieves the calculation parameter struct
cP = retParaStruct(pData.cP);

% removes the dependencies if not separating by region
if (~cP.useRegion)
    % resets the x-dependecy flag
    for i = 1:length(pData.oP.yVar)
        pData.oP.yVar(i).xDep = [];
    end
    
    % resets the x-dependecy flag
    for i = 1:length(pData.oP.xVar)
        pData.oP.xVar(i).Type = 'Other';
    end    

    % resets the statistic fields (from multi-comp to comp)
    isS = find(~cellfun(@isempty,field2cell(pData.oP.yVar,'Stats')));
    for i = isS(:)'
        pData.oP.yVar(i).Stats = {'Comp'};
    end
end
