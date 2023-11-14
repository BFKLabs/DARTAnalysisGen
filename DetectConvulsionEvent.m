% --- detects the convulsion events from the vertical climbing assay
function pData = DetectConvulsionEvent(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);

% sets the function name/type
pData.Name = 'Convulsion Event Detection';
pData.Type = {'Pop','Multi'};
pData.fType = [2 1 1 2];
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
    nExp = length(snTot);
    [pData.hasSP,pData.hasSR,pData.hasRC] = deal(true,nExp==1,false);
    pData.sP = initSpecialPara(snTotL,pData,[],1); 
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
rI.Shape = '1D';
rI.Stim = 'None';
rI.Spec = 'None';

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
hasExD = any(arrayfun(@(x)(~isempty(getExtnData(x))),snTot));
nPara = 5 + hasExD;
cP = setParaFields(nPara);

% sets the tab list names
a = {'1 - Thresholds','2 - Smoothing'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'Number',3,...
    'vEventT','Total Speed Threshold (mm/s)',[0.01 1000 false]);
cP(2) = setParaFields(a{1},'Number',2,...
    'vEventH','Horizontal Speed Threshold (mm/s)',[0.00 1000 false]);        
cP(3) = setParaFields(a{1},'Number',2.5,...
    'yLim','Vertical Height Threshold (mm)',[0.01 1000 false]);
cP(4) = setParaFields(a{1},'Number',0.5,...
    'tEvent','Minimum Duration Threshold (s)',[0.01 1000 false]);        
cP(5) = setParaFields(a{2},'Number',2,...
    'nSm','Smoothing Window Size',[1,20,true]);        
        
% sets the tool-tip strings
cP(1).TTstr = 'Minimum average total speed required for an event to occur';
cP(2).TTstr = 'Minimum average horizontal speed required for an event to occur';
cP(3).TTstr = 'Maximum vertical height for detecting convulsion events';
cP(4).TTstr = 'Minimum duration of convulsion activity';
cP(5).TTstr = 'Smoothing half-window size';

% adds the add external data (if possible)
if hasExD
    cP(end) = setParaFields(a{1},'Boolean',1,'useExD','Use External Data');
    cP(end).TTstr = 'Uses the external data to set analysis duration';
end

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialisations
mltExp = length(snTot) > 1;
nPara = 7 + ~mltExp*3;
pMet = {'Event/Trace Overlay','Individual K-M Curve','Total K-M Curve',...
         'Individual Event Duration','Total Event Duration','Event Count'};
pType = {'Bar Graph','Boxplot'};
pMetT = {'Height','Speed','Horizontal Speed'};

% initialises the parameter struct
pP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General','2 - Bar/Boxplot Parameters','3 - Trace Parameters'};
if mltExp; pMet = pMet(2:end); end

% sets the trace/discrete plot indices
nMet = length(pMet);
[indP,indT] = deal(nMet+(-2:0),1:(nMet-3));

% sets up the general and bar/boxplot plotting parameter fields
pP(1) = setParaFields(a{1},'List',{1,pMet},'pMet','Plot Metric');
pP(2) = setParaFields(a{1},'Boolean',0,'plotGrid','Show Plot Axis Gridlines');
pP(3) = setParaFields(a{2},'List',{1,pType},'pType','Plot Type',[],{1,indP});
pP(4) = setParaFields(a{2},'Number',0.75,'pW',...
                    'Bar Plot Relative Width',[0 1 false],{{1,3},{indP,1}});
pP(5) = setParaFields(a{2},'Boolean',1,'plotErr',...
                    'Show Errorbars/Outlier',[],{1,indP});

% sets the trace type parameter (multi-expt only)
if ~mltExp
    pP(6) = setParaFields(a{3},'List',{1,pMetT},'pMetT','Trace Type',[],{1,1});
end

% sets the other trace parameter fields
iiP = 6 + ~mltExp;
pP(iiP) = setParaFields(a{3},'Number',1,'lWid',...
                'Plot Line Width',[0.1 10 false],{1,indT});

% adds in the trace checkbox parameters            
if ~mltExp
    iiP = 9;
    pP(8) = setParaFields(a{3},'Boolean',0,'useFrm','Plot Frame Index',[],{1,1});            
    pP(9) = setParaFields(a{3},'Boolean',0,'pltThresh','Plot Threshold Limits',[],{1,1});            
end
            
% sets the other trace parameter fields            
pP(iiP+1) = setParaFields(a{3},'Boolean',0,'normSig',...
                'Plot Normalised Signal',[],{1,1+(~mltExp)});            

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
Stats = [];
[Type,Type2,Type3,Type4] = deal(1,[1,5],4,5);
oP = setupOutputParaStruct(snTot);

% sets the independent variable fields
oP = addXVarField(oP,'Event #','indE','Group');
oP = addXVarField(oP,'Time (s)','xCDF','Other');
oP = addXVarField(oP,'Time (s)','xCDFT','Other');

% sets the dependent variable fields
oP = addYVarField(oP,'Count','nEvent',Stats,Type,[],1);
oP = addYVarField(oP,'Total Duration','tEventT',Stats,Type,[],1);
oP = addYVarField(oP,'Event Start','tStart',Stats,Type2,{'indE'},1);
oP = addYVarField(oP,'Event Duration','tEvent',Stats,Type2,{'indE'},1);
oP = addYVarField(oP,'Avg. Speed','vEvent',Stats,Type2,{'indE'},1);
oP = addYVarField(oP,'Count (Indiv)','yCDF',[],Type3,{'xCDF'});
oP = addYVarField(oP,'Count (Total)','yCDFT',[],Type3,{'xCDFT'});

% single experiment dependent variable fields
if length(snTot) == 1
    oP = addXVarField(oP,'Time','T','Time');
    oP = addYVarField(oP,'Height','Y',[],Type4,{'T'},1);
    oP = addYVarField(oP,'Speed','V',[],Type4,{'T'},1);
    oP = addYVarField(oP,'Horz Speed','Vh',[],Type4,{'T'},1);
end

% --- sets the data cursor update function
function dTxt = dataCursorFunc(~,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
pP = retParaStruct(dcObj.pData.pP);
sP = retParaStruct(dcObj.pData.sP);

% sets the common class fields
dcObj.useGrpHdr = false;
dcObj.grpName = dcObj.pData.appName(sP.Sub.isPlot);

% sets the metric specific fields
switch pP.pMet
    case 'Event/Trace Overlay'
        % case is the trace/event overlay
        
        % field retrieval
        uD = evnt.Target.UserData;        
        pPos = dcObj.evnt.Position;                
        
        % sets the plot type based on the selection
        if strcmp(get(evnt.Target,'Marker'),'none')
            dcObj.pType = 'Trace';            
        else
            % sets the plot type
            dcObj.pType = 'Marker';
            
            % determines the index of the selected marker
            hP = evnt.Target.Parent;               
            mType = evnt.Target.Marker;
            hM = findall(hP,'UserData',uD,'Marker',mType);            
            dcObj.mIndex = find(pdist2([hM.XData(:),hM.YData(:)],pPos)==0);
        end
        
        % retrieves the plot data struct
        p = dcObj.plotD{1}(sP.pInd);    
        
        % sets the selected frame index (based on the x-axis type)
        if pP.useFrm
            % case is using the frame index
            iiT = pPos(1);
            dcObj.xName = 'Frame Index';
        else
            % case is using the experiment time
            iiT = p.T{1} == pPos(1);        
            dcObj.xName = 'Time';
            dcObj.tUnits = 'sec';            
        end
        
        % resets the plot metric value
        switch pP.pMetT
            case 'Height'                
                pPos(2) = p.Y{uD}(iiT);
                dcObj.yUnits = 'mm';
                
            case 'Speed'
                pPos(2) = p.V{uD}(iiT);
                dcObj.yUnits = 'mm/s';                
                
            case 'Horizontal Speed'
                pPos(2) = p.Vh{uD}(iiT);
                dcObj.yUnits = 'mm/s';                                
        end
        
        % rescales the y-value
        dcObj.evnt.Position = pPos;
        
        % sets the other properties        
        dcObj.yName = pP.pMetT;        
        [dcObj.useGrpHdr,dcObj.combFig] = deal(true,false);        
        
    case 'Individual K-M Curve'
        % case is the individual K/M curve
        dcObj.pType = 'Trace';
        dcObj.yName = 'Event Count';      
        dcObj.xName = 'Time';
        dcObj.tUnits = 'sec';  
        [dcObj.useGrpHdr,dcObj.combFig] = deal(true);
        
    case 'Total K-M Curve'        
        % case is the total K/M curve
        dcObj.pType = 'Trace';
        dcObj.yName = 'Percentage';
        dcObj.xName = 'Time';
        dcObj.tUnits = 's';
        [dcObj.useGrpHdr,dcObj.combFig] = deal(true);
        
    case 'Event Count'        
        % case is the other plot metrics
        dcObj.pType = pP.pType;
        dcObj.yName = 'Event Count';
        dcObj.xName = 'Group Name';
        dcObj.xGrp = dcObj.grpName;
        
    otherwise
        % case is the other plot metrics
        dcObj.pType = pP.pType;
        dcObj.yName = 'Time';
        dcObj.yUnits = 's';
        dcObj.xName = 'Group Name';
        dcObj.xGrp = dcObj.grpName;        
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
    cP = retParaStruct(pData.cP,gPara);
end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,'T',[],'Y',[],'V',[],'Vh',[],...
                'tEvent',[],'tEvent_mn',[],'tEvent_sem',[],...                
                'tEventT',[],'tEventT_mn',[],'tEventT_sem',[],...                            
                'vEvent',[],'vEvent_mn',[],'vEvent_sem',[],...
                'nEvent',[],'nEvent_mn',[],'nEvent_sem',[],...                            
                'iEvent',[],'tStart',[],'indE',[],'tLim',[],...
                'xCDF',[],'yCDF',[],'xCDFT',[],'yCDFT',[],'xScl',[]);
                        
% creates the waitbar figure
wStr = {'Overall Progress','Event Detection Calculations'};
wOfs = (nExp > 1); wStr = wStr((2-wOfs):end);

% creates the waitbar figure
h = ProgBar(wStr,'Convulsion Event Detection Calculations'); 

% determines the frame count for each video
nFrm0 = arrayfun(@(y)(cellfun(@(x)(length(x)),y.T)),snTot,'un',0);

% --------------------------- %
% --- METRIC CALCULATIONS --- %
% --------------------------- %

% sets up the filter
hF = [ones(cP.nSm+1,1);zeros(cP.nSm,1)];

% retrieves the analysis time limits for each expt
tLim = cell2mat(arrayfun(@(x)(getTimeLimits(x,cP)),snTot(:),'un',0));

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
    
    % sets the acceptance flags
    T0 = cell2mat(snTot(i).T(:));
    flyok = snTot(i).iMov.flyok; 
    
    % determines the feasible data points for the experiment
    iiT = (T0 >= tLim(i,1)) & (T0 <= tLim(i,2));
    
    % calculates the time difference vector
    T = T0(iiT) - tLim(i,1);
    dT0 = imfilter([diff(T0(1:2));diff(T)],hF);
    
    % calculates the convulsion events
    for j = 1:nApp
        % updates the waitbar figure
        wStrNw = sprintf(['Calculating Event Metrics (Region ',...
                          '%i of %i)'],j,nApp);
        if h.Update(1+wOfs,wStrNw,0.50*(1+(j/(nApp+1))))
            % if the user cancelled, then exit the function
            [plotD,ok] = deal([],false);
            return            
        end
        
        % only calculate if values exist...
        if ~isempty(snTot(i).Px{j})          
            % memory allocation
            if isempty(plotD(j).iEvent)
                nFlyMx = size(plotD(j).tStart,2);
                plotD(j).iEvent = cell(1,nFlyMx,nExp);
                plotD(j).indE = cell(nExp,1);
                plotD(j).tLim = tLim;
            end
            
            % calculates the inter-frame speeds
            [V,Vh,X] = calcInterFrameSpeed(snTot(i),j,dT0,hF);
            H = calcPosHeights(X);
            
            % sets the x-coordinates (single expt only)            
            if nExp == 1                  
                plotD(j).T = {T};
                plotD(j).Y(1,flyok{j},i) = num2cell(H(iiT,:),1); 
                plotD(j).V(1,flyok{j},i) = num2cell(V(iiT,:),1);
                plotD(j).Vh(1,flyok{j},i) = num2cell(Vh(iiT,:),1);
            else
                plotD(j).T{i} = T;
            end            
            
            % determines the frames which meet the speed/height tolerances
            [VT,VH,YT] = deal(num2cell(V,1),num2cell(Vh,1),num2cell(H,1));
            iGrp = cellfun(@(x,y,z)(getGroupIndex((x >= cP.vEventT) & ...
                    (y <= cP.yLim) & (z >= cP.vEventH))),VT,YT,VH,'un',0);
               
            % reduces the events to those that meet the duration tolerance
            jGrp = cellfun(@(x)(x(cellfun(@(z)...
                (diff(T0(z([1,end])))),x) >= cP.tEvent)),iGrp,'un',0);

            % sets the start, duration and avg speed of each event
            [tS,tE,vE,iE] = getEventProps(T0,VT,jGrp,nFrm0{i},tLim(i,:));  
            plotD(j).iEvent(1,flyok{j},i) = iE;            
            plotD(j).tStart(1,flyok{j},i) = tS;
            plotD(j).tEvent(1,flyok{j},i) = tE;
            plotD(j).vEvent(1,flyok{j},i) = vE;
            
            % sets the event count and total duration values
            plotD(j).nEvent(1,flyok{j},i) = ...
                            cellfun(@(x)(size(x,1)),iE,'un',0);
            plotD(j).tEventT(1,flyok{j},i) = ...
                            cellfun(@(x)(sum(x,'omitnan')),tE,'un',0);
        end        
    end
    
    % sets the event index array
    nMaxE = max(arrayfun(@(x)(length(x.tStart{1,1,i})),plotD));    
    for j = 1:nApp
        % calculates the array size expansion
        dN = nMaxE - length(plotD(j).tStart{1,1,i});
        if dN > 0
            for k = 1:size(plotD(j).tStart,2)
                % pads the start time arrays
                indP = [k,i];
                plotD(j).tStart = expandDataArray(plotD(j).tStart,indP,dN);
                plotD(j).tEvent = expandDataArray(plotD(j).tEvent,indP,dN);
                plotD(j).vEvent = expandDataArray(plotD(j).vEvent,indP,dN);
            end
        end
        
        % sets the 
        plotD(j).indE{i} = arrayfun(@num2str,(1:nMaxE)','un',0);
        
        % calculates the scale factors
        if nExp == 1
            plotD(j).xScl = zeros(1,3);
            plotD(j).xScl(1) = calcScaleFactor('Height',snTot(i),j);            
            plotD(j).xScl(2) = calcScaleFactor('Speed',plotD(j).V);        
            plotD(j).xScl(3) = calcScaleFactor('Speed',plotD(j).Vh);
        end
    end
    
    % updates the waitbar figure
    h.Update(1+wOfs,'Event Detection Complete!',1);    
end

% % ------------------------- %
% % --- METRIC CONVERSION --- %
% % ------------------------- %
% 
% % determines the maximum event count
% nMax = max(arr2vec(cellfun(@(x)(max(cellfun('length',x))),tS)));
% indES = arrayfun(@(x)(sprintf('Event #%i',x)),(1:nMax),'un',0);
% 
% % metric conversion
% tS = convertVarSizedArrays(tS,nMax);
% tE = convertVarSizedArrays(tE,nMax);
% vE = convertVarSizedArrays(vE,nMax);
% 
% Atmp = NaN(1,nMax);
% for j = 1:nApp
%     % memory allocation 
%     plotD(j).tStart(:) = {Atmp};    
%     plotD(j).tEvent(:) = {Atmp};
%     plotD(j).vEvent(:) = {Atmp};
%     
%     % sets the index string array
%     plotD(j).indE = indES;
%     
%     % sets the raw data values
%     for i = 1:nExp
%         plotD(j) = setRawDataValues(...
%                         plotD(j),snTot(i),tS{i,j},[],'tStart',i,j,2);        
%         plotD(j) = setRawDataValues(...
%                         plotD(j),snTot(i),tE{i,j},[],'tEvent',i,j,2);
%         plotD(j) = setRawDataValues(...
%                         plotD(j),snTot(i),vE{i,j},[],'vEvent',i,j,2);                            
%     end
% end

% ------------------------------- %
% --- STATISTICS CALCULATIONS --- %
% ------------------------------- %

% pre-calculations
xCDF = 0:max(diff(tLim,[],2));

%
tEventT0 = field2cell(plotD,'tEventT');
tEventT = cellfun(@(x)(reduceCellArray(x)),tEventT0,'un',0);
tMaxT = max(cellfun(@(x)(max(cell2mat(arr2vec(x)))),tEventT));
xCDFT = 0:0.01:tMaxT;

% calculates the metrics for all apparatus
for j = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf(['Calculating Event Statistics (Region ',...
                      '%i of %i)'],j,nApp);
    if h.Update(1+wOfs,wStrNw,(j+1)/(2+nApp))
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    end  
    
    % sets up the individial event cdf    
    plotD(j).xCDF = xCDF;    
    Ynn = rmmissing(cell2mat(arr2vec(plotD(j).tStart)));
    if isempty(Ynn)
        plotD(j).yCDF = zeros(size(xCDF))';
    else        
        plotD(j).yCDF = [0;cumsum(histcounts(Ynn,xCDF))'];
    end
             
    % sets up the total duration cdf
    plotD(j).xCDFT = xCDFT;
    tE = rmmissing(arr2vec(cell2mat(tEventT{j})));    
    if isempty(tE)
        plotD(j).yCDFT = zeros(size(xCDFT))';
    else
        yCDFT = [0;cumsum(histcounts(tE(tE > 0),xCDFT))'];
        plotD(j).yCDFT = yCDFT/yCDFT(end);
    end
    
    % calculates the metric statistics for each of the values
    plotD(j) = calcMetricStats(plotD(j),'tEvent',3);
    plotD(j) = calcMetricStats(plotD(j),'vEvent',3);
    plotD(j) = calcMetricStats(plotD(j),'tEventT');
    plotD(j) = calcMetricStats(plotD(j),'nEvent');
end

% closes the waitbar figure
if ~h.Update(1+wOfs,'Event Detection Calculations Complete!',1)
    h.closeProgBar();
end

% --- expands the data arrays
function Y = expandDataArray(Y,iP,dN)

if isempty(Y{1,iP(1),iP(2)})
    Y{1,iP(1),iP(2)} = NaN(dN,1);
else
    Y{1,iP(1),iP(2)} = padarray(Y{1,iP(1),iP(2)},[dN,0],NaN,'post');
end

% ----------------------------------------------------------------------- %
% ---                        PLOTTING FUNCTION                        --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,varargout] = plotFunc(snTot,pData,plotD,ind)

% initialisations
varargout = {};

% retrieves the plotting paraeter struct
pP = retParaStruct(pData.pP);
cP = retParaStruct(pData.cP);
sP = retParaStruct(pData.sP);
pF = pData.pF;

% updates the data cursor
pData.dcFunc = @dataCursorFunc;    
pMetT = {'Height','Speed','Horizontal Speed'};

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the plotting indices and subplot indices
ind = find(sP.Sub.isPlot);
nApp = length(ind); if (nApp == 0); return; end

% sets the plot data struct
if strcmp(pP.pMet,'Event/Trace Overlay')
    p = plotD{1}(sP.pInd);
else
    p = plotD{1}(ind);
end

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% sets the default title/label strings
pF.Title(1).String = pP.pMet;
pF.Title(1).String = pP.pMet;
pF.yLabel(1).String = 'Time (s)';

% sets the other font properties
pF.Axis.Font.FontSize = 16;

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% sets up the figure
hFig = get(hP,'Parent');
hAx = createSubPlotAxes(hP);
axis(hAx,'on'); hold on; 

% sets the plot values/axis labels based on the metric type
switch pP.pMet
    case 'Event/Trace Overlay'
        % case is the trace/event overlay
        
        % initialisations
        xiF = find(snTot.iMov.flyok{sP.pInd});
        nFly = xiF(end);                
        
        % sets the time vector
        if pP.useFrm
            xLblStr = 'Frame Index';
            TT = 1:length(p.T{1});
        else
            xLblStr = 'Time (sec)';
            TT = p.T{1};
        end
        
        % plot object properties
        iPltM = find(strcmp(pMetT,pP.pMetT));
        [yLim,xLim] = deal([0,nFly],[0,TT(end)]);
        yTickLbl = arrayfun(@num2str,1:nFly,'un',0);
        tStr = sprintf('%s (%s)',pP.pMetT,pData.appName{sP.pInd});        
        xScl = p.xScl(iPltM);        
        
        % plots the traces        
        switch pP.pMetT
            case 'Height'
                yThresh = cP.vEventT;
                Yplt = cellfun(@(x,y)(scaleSig(x/xScl,y-1)),...
                                    p.Y(xiF),num2cell(xiF)','un',0);                                

            case 'Speed'
                yThresh = cP.yLim;                
                Yplt = cellfun(@(x,y)(scaleSig(x/xScl,y-1)),...
                                    p.V(xiF),num2cell(xiF)','un',0);
                                
            case 'Horizontal Speed'
                yThresh = NaN;
                Yplt = cellfun(@(x,y)(scaleSig(x/xScl,y-1)),...
                                    p.Vh(xiF),num2cell(xiF)','un',0);                                
        end
        
        % plots the tracks
        xiY = num2cell(1:length(Yplt));
        cellfun(@(x,y)(plot(TT,x,'b','UserData',y)),Yplt,xiY)                
        
        % sets the event/separation lines
        for i = 1:nFly
            % plots the event markers (if available)
            if ~isempty(p.tStart{i})
                % determines the time points indices
                Ts = arr2vec(p.tStart{i});
                [~,indP,iC] = intersect(p.T{1},Ts);
                
                % creates the plot object
                if pP.useFrm
                    plot(indP,Yplt{i}(indP),'r.','tag','hEvent',...
                                'UserData',i,'MarkerSize',20);                    
                else
                    plot(Ts(iC),Yplt{i}(indP),'r.','tag','hEvent',...
                                'UserData',i,'MarkerSize',20);
                end
                
                % plots the threshold marker
                if pP.pltThresh
                    yPltL = scaleSig(yThresh/xScl,i-1)*[1,1];
                    plot(TT([1,end]),yPltL,'k--');
                end
            end
            
            % plots the separation line
            if i < nFly
                plot(xLim,i*[1,1],'k')
            end
        end        
            
        % sets the text label
        hText = imtext(0,0,{''},'right');
        set(hText,'parent',hAx,'visible','off','FontSize',8,...
                'FontWeight','bold','EdgeColor','k','LineWidth',1,...
                'BackgroundColor','y','tag','hText');        
        
        % sets the axis properties
        set(hAx,'ylim',yLim,'xLim',xLim,...
                'yTick',(1:nFly)-0.5,'yTickLabel',yTickLbl);        
            
        % sets the title/label strings
        pF.Title(1).String = tStr;
        pF.yLabel(1).String = 'Fly Index';
        pF.xLabel(1).String = xLblStr;
        
        % flips the vertical axis
        axis(hAx,'ij')
        
    case 'Individual K-M Curve'
        % case is the individual Kaplan-Meier curve

        % memory allocation
        hPlt = zeros(nApp,1);
        xLim = [0,max(diff(p(1).tLim,[],2))];
        
        % creates the K-M curves
        for i = 1:nApp
            [x,f] = deal(p(i).xCDF,p(i).yCDF);
            if pP.normSig
                hPlt(i) = stairs(hAx,x,f/f(end),'LineWidth',pP.lWid,...
                    'UserData',i);
            else
                hPlt(i) = stairs(hAx,x,f,'LineWidth',pP.lWid,'UserData',i);
            end
        end

        % creates the legend
        pF.Legend.String = pData.appName(ind);
        hLg = createLegendObj(hPlt,pF.Legend);
        set(hLg,'Location','NorthWest');           
        
        % sets the y-axis label string
        if pP.normSig
            pF.yLabel(1).String = 'Proportion';
        else
            pF.yLabel(1).String = 'Event Count';
        end
        
        % sets the axis labels and other properties
        pF.xLabel(1).String = 'Time (sec)';
        set(hAx,'xlim',xLim)
        
    case 'Total K-M Curve'
        % case is the total duration Kaplan-Meier curve
        
        % memory allocation
        hPlt = zeros(nApp,1);            
        
        % creates the survival curves for each genotype group
        for i = 1:nApp
            [x,f] = deal(p(i).xCDFT,p(i).yCDFT);
            hPlt(i) = stairs(hAx,x,100*(1-f),'LineWidth',pP.lWid,...
                'UserData',i);
        end        
        
        % creates the legend
        pF.Legend.String = pData.appName(ind);
        hLg = createLegendObj(hPlt,pF.Legend);
        set(hLg,'Location','NorthEast'); 
        
        % sets the title/label strings
        pF.yLabel(1).String = 'Percentage';
        pF.xLabel(1).String = 'Event Duration (s)';
        set(hAx,'xlim',[cP.tEvent,max(get(hAx,'xlim'))],'ylim',[0,100])
        
    otherwise
        % individual/total event duration
        
        % initialisations
        xi = 1:nApp;                
        
        % metric specific initialisations
        switch pP.pMet
            case 'Individual Event Duration'
                pStr = 'tEvent';
            case 'Total Event Duration'
                pStr = 'tEventT';
            case 'Event Count'
                pStr = 'nEvent';
                pF.yLabel(1).String = 'Convulsion Events';
        end        
        
        % creates the bar/boxplot metrics
        plotBarBoxMetrics(hAx,xi,p,pStr,pP,[],'b',1);       
        set(hAx,'XTickLabel',pData.appName(ind))
        
end
    
% formats the plot axis
formatPlotAxis(hAx,pF,1);

% updates the main window motion function
hToolDC = findall(hFig,'tag','menuDataCursor');
if strcmp(get(hToolDC,'state'),'on')
    set(hFig,'WindowButtonMotionFcn',pData.dcFunc);
end

% sets the other axis properties
if pP.plotGrid; grid on; end
set(hAx,'Box','on')

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% finish me?!
a = 1;
   
% ----------------------------------------------------------------------- %
% ---                         OTHER FUNCTIONS                         --- %
% ----------------------------------------------------------------------- %

% --- calculates the interframe speed
function [V,Vh,X] = calcInterFrameSpeed(snTot,iApp,dT,hF)

% calculates the squared x-displacements
X = snTot.Px{iApp};
dP2 = diff([X(1,:);X],[],1).^2;

% appends the y-displacement squared values (if available)
if ~isempty(snTot.Py)
    Y = snTot.Py{iApp};
    dY = diff([Y(1,:);Y],[],1);
    Vh = imfilter(abs(dY),hF)./dT;    
    
    dP2 = dP2 + dY.^2;
else
    Vh = NaN(size(dT));
end

% calculates the speed over all frames
V = imfilter(sqrt(dP2),hF)./dT;

% --- calculates the height position of the flies
function Y = calcPosHeights(Px)

% parameters
pTile = 99;

% calculates the height of the fly
Y = max(0,repmat(prctile(Px,pTile,1),size(Px,1),1) - Px);

% --- calculates the event properties
function [tS,tE,vE,iE] = getEventProps(T,VT,jGrp,nFrm,tLim)

% memory allocation
nGrp = cellfun('length',jGrp);
iE = cell(1,length(jGrp));

if max(nGrp) == 0
    [tS,tE,vE] = deal(repmat({NaN},1,length(jGrp)));
    return
else
    [tS,tE,vE] = deal(repmat({NaN(max(nGrp),1)},1,length(jGrp)));
    nFrmT = [0;cumsum(nFrm(1:end-1))];
end

% determines the event properties for each fly
for i = find(~cellfun('isempty',jGrp(:)')) 
    % sets the indices of the event start points
    iS = cellfun(@(x)(x(1)),jGrp{i});
    T0 = T(iS);
    
    % determines the time events that fall within the time limits
    jj = (T0 >= tLim(1)) & (T0 <= tLim(2));        
    if any(jj)
        % if so, set the event properties
        ii = 1:sum(jj);    
        tS{i}(ii) = T0(jj) - tLim(1);   
        tE{i}(ii) = cellfun(@(x)(diff(T(x([1,end])))),jGrp{i}(jj));
        vE{i}(ii) = cellfun(@(x)(mean(VT{i}(x))),jGrp{i}(jj));

        % sets the event video/frame indices
        iVid = arrayfun(@(x)(find(x>=nFrmT,1,'last')),iS(jj));
        iFrm = arrayfun(@(x,iv)(x-sum(nFrmT(1:iv))),iS(jj),iVid);
        iE{i} = [iVid,iFrm];
    end
end

% reduces down the arrays to only include actual events
xiF = 1:max(1,max(cellfun(@(x)(size(x,1)),iE)));
tS = cellfun(@(x)(x(xiF)),tS,'un',0);
tE = cellfun(@(x)(x(xiF)),tE,'un',0);
vE = cellfun(@(x)(x(xiF)),vE,'un',0);

% --- removes the empty cells from the array, Y
function Y = reduceCellArray(Y)

Y = arr2vec(Y(~cellfun('isempty',Y)));

% --- calculates the scaled signal from the base signal, Y0
function Ys = scaleSig(Y0,yOfs)

% parameters
pW = 0.025;

% calculates the scaled signal
Ys = (yOfs + pW) + (1-2*pW)*(1-Y0);

% --- gets the analysis time limits
function tLim = getTimeLimits(snTot,cP)

% initialisations
pStr = 'TimingHS';
hStr = {'Start (s)'  'Duration (s)'};

% sets the full experiment time limit
T = cell2mat(snTot.T(:));
tLim = arr2vec(T([1,end]))';

% if there is no external data (or not using) then exit the function
if isempty(snTot.exD) || ~cP.useExD; return; end

% determines if there are any matching external data fields
pStrS = cellfun(@(x)(x.pStr),snTot.exD,'un',0);
hStrS = cellfun(@(x)(x.hStr),snTot.exD,'un',0);
iiP = cellfun(@(x,y)(strcmp(pStr,x) && isequal(hStr,y)),pStrS,hStrS);

% if so, then retrieve the relevant data field
if any(iiP)
    Data = cellfun(@str2double,snTot.exD{iiP}.Data);
    if ~any(isnan(Data))
        % if the data is feasible, then set the time limit
        tLim = Data(1) + [0,Data(2)];
    end
end

% --- calculates the metric scale factor (based on type)
function xScl = calcScaleFactor(pMetT,sObj,varargin)

switch pMetT
    case 'Height'
        % case is the height 
        iC = sObj.iMov.iC{varargin{1}};
        xScl = range(iC)*sObj.sgP.sFac;

    case 'Speed'
        % case is the speed
        xScl = 1.05*max(cellfun(@max,sObj(~cellfun('isempty',sObj))));        
end


% --- retrieves the position coordinates (smooths is necessary)
function Y = getPosCoords(Y,cP)

if cP.useSm
    nSm = 2*cP.nSm + 1;
    Y = cell2mat(cellfun(@(x)(smooth(x,nSm)),num2cell(Y,1),'un',0));
end

% --- calculates the central difference calculations
function dZ = calcBackwardDiff(Z)

% memory allocation
dZ = zeros(size(Z));

% calculates the backward differences
dZ(2:end,:) = Z(2:end,:) - Z(1:(end-1),:);

% --- calculates the central difference calculations
function dZ = calcCentralDiff(Z)

% memory allocation
dZ = zeros(size(Z));

% calculates the central derivatives (end points are forward difference)
dZ(1,:) = diff(Z(1:2));
dZ(2:end-1,:) = Z(3:end,:) - Z(1:(end-2),:);
dZ(end,:) = diff(Z((end-1):end,:));