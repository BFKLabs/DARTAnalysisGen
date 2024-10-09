% --- plot the x/y-location of the fly over the duration of the experiment 
%     (2D analysis only)
function pData = FlyPositionTrace2D(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);

% sets the function name/
pData.Name = '2D Fly Position Traces';
pData.Type = {'Indiv','Multi'};
pData.fType = [2 1 1 3];
pData.rI = initFuncReqInfo(pData);

% initialises the other fields  (if input argument provided)
if (nargin == 1)
    % parameter data struct initialisation
    snTotL = snTot(1);
    pData.cP = initCalcPara(snTot);
    pData.pP = initPlotPara(snTot);
    pData.oP = initOutputPara(snTot);
    pData.pF = initPlotFormat(snTotL);
    
    % sets up the special parameter flags
    [pData.hasRC,pData.hasRS] = deal(false);        
    
    if detMltTrkStatus(snTotL.iMov)
        % case is multi-tracking
        
        % sets the region name/count
        pInfo = snTotL.iMov.pInfo;
        xiR = (1:pInfo.nRow*pInfo.nCol)';
        pData.appName = arrayfun(@(x)(sprintf('Region #%i',x)),xiR,'un',0);
        pData.nApp = length(pData.appName);           
        
        % special parameter fields/data struct
        [pData.hasSR,pData.useReg] = deal(true);
        pData.sP = initSpecialPara(snTotL,pData,[],1);        
               
    else
        % case is single tracking
        
        % sets the apparatus name/count
        pData.appName = snTotL.iMov.pInfo.gName;
        pData.nApp = length(pData.appName);        
        
        % special parameter fields/data struct
        pData.hasSP = true;
        pData.sP = initSpecialPara(snTotL,pData);
    end
    
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
rI.Shape = '2D';
rI.Stim = 'None';
rI.Spec = 'None';
        
% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 3;
cP = setParaFields(nPara);

% sets the calculation parameter fields into the data struct
cP(1) = setParaFields([],'Boolean',0,'useAll','Analyse Entire Experiment');
cP(2) = setParaFields([],'Number',0,'T0','Start Time (min)',[0 inf true],{1,1});
cP(3) = setParaFields([],'Number',10,'Tdur','Analysis Duration (min)',[1 inf true],{1,1});

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 7;
pP = setParaFields(nPara);

% sets the count string (based on tracking type)
if detMltTrkStatus(snTot.iMov)
    cStr = 'Plot Fly Count';
else
    cStr = 'Plot Row Count';
end

% parameter tab strings
a = {'1 - General','2 - Trace Markers'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},'Number',10,'nRow',cStr,[1 100 true]);
pP(2) = setParaFields(a{1},'Number',1,'lWid','Plot Line Width',[0.1 10 false]);
pP(3) = setParaFields(a{1},'Boolean',1,'pltAvg','Plot Average Traces');
pP(4) = setParaFields(a{1},'Boolean',0,'randPerm','Traces Selected Randomly',[],{3,1});
pP(5) = setParaFields(a{2},'Boolean',1,'showStart','Show Start Marker');
pP(6) = setParaFields(a{2},'Boolean',1,'showFinish','Show Finish Marker');
pP(7) = setParaFields(a{1},'Boolean',1,'showTitles','Show Genotype Titles');

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

% checks to see if the solution struct has the sub-region data struct
ok = checkFuncPara({'HasSubRegionStruct'},cP,snTot);
if ~ok; plotD = []; return; end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensions
isMT = detMltTrkStatus(snTot.iMov);
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                'T',T,'X',[],'Y',[],'Davg',[],'iFrm',[],'iFly',[]);

% other memory allocations
[Px,Py] = deal(cell(1,nApp));

% ------------------------------------------------ %
% --- TIME BIN THRESHOLDING & POSITION SETTING --- %
% ------------------------------------------------ %

% creates the waitbar figure
wStr = {'Overall Progress'};
h = ProgBar(wStr,'Movement Trace Calculations');

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
    
    % sets the relevant time points and apparatus indices for this expt
    if cP.useAll
        % uses all the time points
        ii = 1:length(T{i});
    else
        % use only the points from the start to the duration end
        ii = (T{i} >= 60*cP.T0) & (T{i} <= 60*(cP.T0 + cP.Tdur));
    end
    
    % calculates the video frame rate and experiment apparatus indices
    iApp = find(~cellfun('isempty',snTot(i).iMov.flyok));
    
    % sets the relevant x/y-locations for the current experiment  
    [dPx,dPy,R] = get2DCoordsBG(snTot(i),iApp,ii);
    if (length(iApp) == 1); [dPx,dPy,R] = deal({dPx},{dPy},{R}); end
    
    % determines which trace belongs to which sub-region. from this, scale
    % the data values
    for j = 1:length(iApp)
        % calculates the normalised x/y-locations        
        switch snTot(i).iMov.pInfo.mShape
            case 'Circle'
                % case is circular regions
                RR = repmat(R{j}(:)',size(dPx{j},1),1);
                Px{iApp(j)} = [Px{iApp(j)},num2cell(0.5+dPx{j}./(2*RR),1)];
                Py{iApp(j)} = [Py{iApp(j)},num2cell(0.5+dPy{j}./(2*RR),1)];
                
            case 'Rectangle'
                % case is rectangular regions
                [W,H] = deal(R{j}(1,:),R{j}(2,:));
                Px{iApp(j)} = [Px{iApp(j)},num2cell(0.5+dPx{j}./W,1)];
                Py{iApp(j)} = [Py{iApp(j)},num2cell(0.5+dPy{j}./H,1)];
        end
    end
end

% sets trace values for each apparatus
for i = 1:nApp
    if isMT
        % case is multi-tracking    
        
        % determines the indices of the flies from each region
        cIDR = snTot.cID{i};
        [~,~,iC] = unique(cIDR(:,1),'stable');
        indC = arrayfun(@(x)(find(iC == x)),1:max(iC),'un',0);

        % separates the flies by region
        for k = 1:length(indC)
            % sets the x/y-coordinates
            j = cIDR(indC{k}(1),1);
            plotD(j).X = [plotD(j).X,Px{i}(indC{k})];
            plotD(j).Y = [plotD(j).Y,Py{i}(indC{k})];           
            
            if i == nApp
                % calculates the total displacement
                plotD(j).Davg = cellfun(@(x,y)(...
                    calcTotalDist(x,y)),plotD(j).X,plotD(j).Y);

                % combines the coordinate arrays
                plotD(j).X = cell2mat(plotD(j).X);
                plotD(j).Y = cell2mat(plotD(j).Y);
                
                % sets the fly index strings
                plotD(j) = setupIndexString(plotD(j));
            end
        end
        
    else
        % case is single tracking
        
        % sets the x/y-position and average speed values
        plotD(i).X = combineNumericCells(Px{i});    
        plotD(i).Y = combineNumericCells(Py{i});    
        plotD(i).Davg = cellfun(@(x,y)(calcTotalDist(x,y)),Px{i},Py{i});
                                        
        % sets the fly indices
        plotD(i) = setupIndexString(plotD(i));
    end
end

% closes the waitbar figure
if ~h.Update(1,'Activity Calculations Complete!',1)
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

% field flags
isMT = detMltTrkStatus(snTot.iMov);

% sets the plotting indices and subplot indices
if isMT
    % case is multi-tracking
    [p,nApp] = deal(plotD{1}(sP.pInd),1);
    
else
    % case is single tracking
    ind = find(sP.Sub.isPlot);
    nApp = length(ind); if (nApp == 0); return; end
    p = plotD{1}(ind);
end

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');
delete(getCurrentAxesProp)

% memory allocation
hAx = zeros(nApp,1);

% sets the plot line colours 
col = num2cell(distinguishable_colors(nApp,'w'),2);

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% resets the title strings
if isMT
    % case is multi-tracking
    pF.Title.String = sprintf('Region %i',sP.pInd);
    
else
    % case is single-tracking
    pF.Title = pF.Title(ind);
end

% retrieves the formatting struct
pF = retFormatStruct(pF,nApp);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% retrieves the X/Y values
dH = 0.075;
[X,Y] = field2cell(p,{'X','Y'});
pP.nRow = min(pP.nRow,max(cellfun(@(x)(size(x,2)),X)));

% loops through all the subplots
for i = 1:nApp
    % creates the subplot axis
    hAx(i) = axes('OuterPosition',calcOuterPos(1,nApp,i),'parent',hP); 
    resetObjPos(hAx(i),'Bottom',-dH,1);
    resetObjPos(hAx(i),'Height',dH,1);
    hold(hAx(i),'on');    
        
    % determines the optimal configuration
    if i == 1 
        [nRow,nCol] = detMostSquareSetup(hAx(i),pP.nRow);
        [delX,delY] = deal(1e-4*nCol,1e-4*nRow);
    end    
        
    % sets the actual plot index and creates the subplot region
    if size(X{i},2) < pP.nRow
        iPlt = 1:size(X{i},2);
    else
        if pP.pltAvg
            % plots the traces that are closest to the average displacement
            [~,a] = sort(abs(p(i).Davg - median(p(i).Davg)),'ascend');
            iPlt = a(1:pP.nRow); iPlt = iPlt(randperm(pP.nRow));
        elseif (pP.randPerm)
            % randomly permute the plot indices
            iPlt = sort(randperm(size(X{i},2),pP.nRow),'ascend');
        else
            % set the first nRow indices
            iPlt = 1:pP.nRow;
        end
    end
    
    % plots the traces
    for j = 1:min(length(iPlt),pP.nRow)
        % sets the new plot values
        [xOfs,yOfs] = deal(mod(j-1,nCol),(nRow-1)-floor((j-1)/nCol));
        [xPltNw,yPltNw] = deal(X{i}(:,iPlt(j))+xOfs,Y{i}(:,iPlt(j))+yOfs);
        
        % plots the 2D traces 
        plot(xPltNw,yPltNw,'color',col{i},'linewidth',pP.lWid);
        
        % plots that start marker (if required)
        if pP.showStart
            i0 = find(~isnan(xPltNw),1,'first');
            plot(xPltNw(i0),yPltNw(i0),'ko','markersize',10,'linewidth',2); 
        end
        
        % plots that finish marker (if required)        
        if pP.showFinish
            i1 = find(~isnan(xPltNw),1,'last');
            plot(xPltNw(i1),yPltNw(i1),'kx','markersize',10,'linewidth',2)
        end
    end
    
    % sets the x/y axis limits
    set(hAx(i),'TickLength',[0 0],'xlim',[0 nCol]+delX*[-1 1],...
               'yLim',[0 nRow]+delY*[-1 1],'xticklabel',[],'yticklabel',[])        
    
    % plots the column seperators
    for j = 0:nRow
        hCol = plot([0 nCol],j*[1 1],'k:');
        if (any(j == [0 nRow])); set(hCol,'linestyle','-','linewidth',1.5); end
    end    

    % plots the row seperators
    for j = 0:nCol
        hCol = plot(j*[1 1],[0 nRow],'k:');
        if (any(j == [0 nCol])); set(hCol,'linestyle','-','linewidth',1.5); end
    end
               
    % formats the plot axis
    formatPlotAxis(hAx(i),pF,i);
end

% optimises the title placements
optTitlePlacement(hAx,'Title')

% ensures the axis is equal is size and hold is on
cellfun(@(x)(axis(x,'equal')),num2cell(hAx))
try
    cellfun(@(x)(set(x,'xcolor','none','ycolor','none')),num2cell(hAx))
catch
    cellfun(@(x)(set(x,'xcolor','w','ycolor','w')),num2cell(hAx))
end

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)



% ----------------------------------------------------------------------- %
% ---                     MISCELLANEOUS FUNCTION                      --- %
% ----------------------------------------------------------------------- %

% --- sets up the fly/frame index strings
function pD = setupIndexString(pD)

% initialisations
[xiF1,xiF2] = deal(1:size(pD.X,2),1:size(pD.X,1));

% sets up the fly/frame indices
pD.iFly = arrayfun(@(x)(sprintf('Fly #%i',x)),xiF1,'un',0);                                             
pD.iFrm = arrayfun(@(x)(sprintf('Frame #%i',x)),xiF2,'un',0);

% --- calculates the total displacement
function D = calcTotalDist(x,y)

D = sum(sqrt(diff(x).^2 + diff(y).^2));