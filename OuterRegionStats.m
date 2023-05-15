% --- calculates the fly activity initiation (number of starts per second 
%     stopped) over the duration of an experiment (short experiment only)
function pData = OuterRegionStats(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Outer Region Statistics';
pData.Type = {'Pop','Multi'};
pData.fType = [2 1 1 3];
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
            'Stim',[],'Spec',[],'SpecFcn',[],'ClassicFcn',false);
        
% sets the struct fields
rI.Scope = setFuncScopeString(pData.Type);
rI.Dur = 'None';
rI.Shape = '2D (Circle)';
rI.Stim = 'None';
rI.Spec = 'None';

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 1;                       
cP = setParaFields(nPara);

% sets the parameter fields
cP(1) = setParaFields([],'Number',3,'rTol','Outer Edge Distance (mm)',[0.5 5 false]);

% sets the tool-tip strings
cP(1).TTstr = 'Radial distance from circle edge encompassing the outer region';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 5;
pP = setParaFields(nPara);

% sets the parameter list strings
pList = {'Outer Edge Crossings','Outer Edge Percentage',...
         'Mean Radial Distance'};
pList2 = {'Boxplot','Bar Graph'};

% parameter tab strings
a = {'1 - General','2 - Metrics'};

% sets the parameter fields
pP(1) = setParaFields(a{1},'List',{1,pList},'pMet','Plot Metric');
pP(2) = setParaFields(a{1},'Boolean',1,'plotGrid','Show Axis Gridlines');
pP(3) = setParaFields(a{1},'Boolean',1,'plotErr','Show Error Bars/Outliers');
pP(4) = setParaFields(a{2},'List',{1,pList2},'pType','Graph Plot Type');
pP(5) = setParaFields(a{2},'Number',0.75,'pW','Bar Plot Relative Width',[0 1 false],{4,2});

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
oP = setupOutputParaStruct(snTot);
[Type1,Type2,Stats] = deal(3,2,{'Comp'});

% sets the dependent output variables
oP = addYVarField(oP,'Outer Edge Crossing','nX',Stats,Type1,[],1);
oP = addYVarField(oP,'Outer Edge Crossing (Mean)','nX_mn',[],Type2);
oP = addYVarField(oP,'Outer Edge Crossing (SEM)','nX_sem',[],Type2);
oP = addYVarField(oP,'Outer Edge Proportion','pOut',Stats,Type1,[],1);
oP = addYVarField(oP,'Outer Edge Proportion (Mean)','pOut_mn',[],Type2);
oP = addYVarField(oP,'Outer Edge Proportion (SEM)','pOut_sem',[],Type2);
oP = addYVarField(oP,'Mean Radial Distance','R',Stats,Type1,[],1);
oP = addYVarField(oP,'Mean Radial Distance (Mean)','R_mn',[],Type2);
oP = addYVarField(oP,'Mean Radial Distance (SEM)','R_sem',[],Type2);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
pP = retParaStruct(dcObj.pData.pP);
sP = retParaStruct(dcObj.pData.sP);

% other initialisations
grpName = dcObj.pData.appName(sP.Sub.isPlot);

% sets the common class fields
dcObj.pType = pP.pType;
dcObj.yName = pP.pMet;
dcObj.grpName = grpName;
[dcObj.useGrpHdr,dcObj.combFig] = deal(false);
[dcObj.xName,dcObj.xGrp] = deal('Group Name',grpName);

% sets up the axes specfic fields
switch pP.pMet
    case 'Outer Edge Crossings'
        % case is the outer edge crossings
        dcObj.yUnits = 'crossings';
        
    case 'Outer Edge Percentage'
        % case is the outer edge percentage
        dcObj.yUnits = '%';        
        
    case 'Mean Radial Distance'
        % case is the mean radial distance
        dcObj.yUnits = 'mm';        
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
cP.movType = 'Absolute Speed';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
[nApp,nExp,ok] = deal(length(snTot(1).iMov.flyok),length(snTot),true);

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                              'nX',[],'nXT',[],'nX_mn',[],'nX_sem',[],...
                              'pOut',[],'pOutT',[],'pOut_mn',[],'pOut_sem',[],...
                              'R',[],'RT',[],'R_mn',[],'R_sem',[]);                                     
                          
% waitbar figure setup
wStr = {'Calculating Relative Locations','Edge Region Detection'};
h = ProgBar(wStr,'Outer Region Statistics');          

% ------------------------------------------------------- %
% --- INTER-FRAME DISTANCE CALCULATION & THRESHOLDING --- %
% ------------------------------------------------------- %

for i = 1:nExp
    % updates the waitbar figure
    wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
    if h.Update(1,wStrNw,i/nExp)
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    else
        h.Update(2,wStr{2},0);
    end
    
    [Texp,sFac] = deal(cell2mat(snTot(i).T),snTot(i).sgP.sFac);
    iApp = find(~cellfun('isempty',snTot(i).iMov.flyok));
    
    for j = 1:length(iApp)
        % updates the waitbar figure
        wStrNw = sprintf('%s (Apparatus %i of %i)',wStr{2},j,length(iApp));
        if h.Update(2,wStrNw,j/length(iApp))
            % if the user cancelled, then exit the function
            [plotD,ok] = deal([],false);
            return            
        end
        
        % retrieves the relative 2D coordinates
        k = iApp(j);
        [dPx,dPy,R] = get2DCoordsBG(snTot(i),k); 

        % determines the edge position flags for each fly (over all frames)
        onEdge = detFlyEdgePos(dPx,dPy,R,cP,sFac);                 
        
        % calculates the outer region statistics for all flies
        rPosM = cellfun(@(x,y,z)(calcOuterRegionStats(Texp,x,y,z,sFac)),...
                num2cell(onEdge,1),num2cell(dPx,1),num2cell(dPy,1),'un',0);
        [nX,pOut,R] = field2cell(cell2mat(rPosM),{'nX','prOut','R'},1);
        
        % sets the values into the plotting data struct
        kk = 1:length(nX);
        plotD(k).nX(1,kk,i) = num2cell(nX,1);
        plotD(k).pOut(1,kk,i) = num2cell(pOut,1);
        plotD(k).R(1,kk,i) = num2cell(R,1);
    end
end

% sets/calculates the raw, mean and SEM proportional activity values
for i = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf('%s (Apparatus %i of %i)',wStr{2},i,nApp);
    if h.Update(2,wStrNw,(2*i-1)/(2*nApp+1))
        [plotD,ok] = deal([],false);
        return
    end    
            
    % calculates the metric statistics
    plotD(i) = calcMetricStats(plotD(i),'nX'); 
    plotD(i) = calcMetricStats(plotD(i),'pOut');
    plotD(i) = calcMetricStats(plotD(i),'R');  
end

% closes the waitbar figure
if ~h.Update(2,'Outer Region Statistic Calculations Complete!',1)
    h.closeProgBar(); 
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

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% retrieves the formatting struct
pF = retFormatStruct(pF,1);
[xi,yL] = deal(1:length(ind),[]);
xLim = xi([1 end])+0.5*[-1 1];    

% sets the title string
pF.Title.String = pP.pMet;

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% sets up the plot axis
hAx = createSubPlotAxes();    
axis(hAx,'on'); hold on; 

% sets the plot variable string and y-axis string
switch (pP.pMet)        
    case ('Outer Edge Crossings')
        [pF.yLabel.String,pStr] = deal('Crossings','nX');        
    case ('Outer Edge Percentage')
        [pF.yLabel.String,pStr,yL] = deal('Percentage','pOut',[0 100]);
    case ('Mean Radial Distance')
        [pF.yLabel.String,pStr] = deal('Distance (mm)','R');
end

% plots the bar/boxplot graph
plotBarBoxMetrics(hAx,xi,p,pStr,pP,yL,'b');

% updates the axis properties
set(hAx,'xticklabel',[],'xlim',xLim,'linewidth',1.5,'UserData',1);

% adds in the gridlines (if checked)
if (pP.plotGrid); grid(hAx,'on'); end

% ------------------------------ %
% --- PLOT AXES REFORMATTING --- %
% ------------------------------ %

% sets the x-axis labels
formatPlotAxis(hAx,pF,1);

% formats and resets the axis positions
resetAxesPos(hAx,1,1); 

% sets the group strings
setGroupString(hAx,pF,xi,snTot(1).iMov.pInfo.gName(ind),30);

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% ----------------------------------------------------------------------- %
% ---                         OTHER FUNCTIONS                         --- %
% ----------------------------------------------------------------------- %

% --- determines the relative positional metrics
function rPosM = calcOuterRegionStats(T,onEdge,X,Y,sFac)

% determines the index groups of the points within the inner/outer regions
rPosM = struct('nX',NaN,'prOut',NaN,'R',NaN);
[inGrp,outGrp] = deal(getGroupIndex(~onEdge),getGroupIndex(onEdge));

% determines the number of outer region crossings and proportion of the
% time where the fly is located within the inner region
rPosM.nX = (length(inGrp) + length(outGrp)) - 1;

% calculates the time when the fly is in the inner/outer regions (for each
% event when this occurs)
tIn = sum(cellfun(@(x)(diff(T(x([1 end])))),inGrp));
tOut = sum(cellfun(@(x)(diff(T(x([1 end])))),outGrp));
rPosM.prOut = 100*tOut/(tIn+tOut);

% calculates the mean radial distance of the fly over the experiment
rPosM.R = mean(sqrt(X.^2 + Y.^2)*sFac,'omitnan');
