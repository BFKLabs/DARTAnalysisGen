% --- calculates the level of fly arousability with respect to delivered 
%     stimuli throughout the day (long experiment only)
function pData = RegionExploration(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);

% sets the function name/type
pData.Name = '2D Region Exploration';
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
rI.Shape = '2D';
rI.Stim = 'None';
rI.Spec = 'None';

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 0;
cP = setParaFields(nPara);

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 1;
pP = setParaFields(nPara);

% sets the plot parameter fields into the data struct
pP(1) = setParaFields([],'Boolean',1,'plotGrid','Show Axis Gridlines');

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
pF = setFormatFields(1);

% initialises the font structs
pF.Title = setFormatFields([],'Exploration Percentage',1);
pF.xLabel = setFormatFields([],'');
pF.yLabel = setFormatFields([],'% Explored');
pF.Axis = setFormatFields([],[]);

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);
[Type1,Type2,Stats] = deal(3,2,{'Comp'});

% sets the dependent output variables
oP = addYVarField(oP,'Exploration %age','E',Stats,Type1,[],1);
oP = addYVarField(oP,'Exploration %age (Mean)','E_mn',[],Type2);
oP = addYVarField(oP,'Exploration %age (SEM)','E_sem',[],Type2);
oP = addYVarField(oP,'Exploration %age (Median)','E_md',[],Type2);
oP = addYVarField(oP,'Exploration %age (Lower Quartile)','E_lq',[],Type2);
oP = addYVarField(oP,'Exploration %age (Upper Quartile)','E_uq',[],Type2);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
sP = retParaStruct(dcObj.pData.sP);

% other initialisations
grpName = dcObj.pData.appName(sP.Sub.isPlot);

% sets the common class fields
dcObj.pType = 'Boxplot';
dcObj.yName = 'Outer Edge Crossings';
dcObj.yUnits = 'crossings';
dcObj.grpName = grpName;
[dcObj.useGrpHdr,dcObj.combFig] = deal(false);
[dcObj.xName,dcObj.xGrp] = deal('Group Name',grpName);

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

% checks to see if the solution struct has the sub-region data struct
snTotL = snTot(1);
ok = checkFuncPara({'HasSubRegionStruct'},cP,snTotL);
if (~ok); plotD = []; return; end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensions
mShape = snTot(1).iMov.pInfo.mShape;
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                              'E',[],'E_mn',[],'E_sem',[],...
                              'E_md',[],'E_lq',[],'E_uq',[]);    

% other memory allocations
DL = [];

% --------------------------------------- %
% --- REGION EXPLORATION CALCULATIONS --- %
% --------------------------------------- %

% creates the waitbar figure
wStr = {'Overall Progress','Group Progress'};
wStr = wStr(1:(1+(nExp > 1)));
h = ProgBar(wStr,'Region Exploration Calculations');

% loops through each of the experiments calculating the velocity values
for i = 1:nExp
    % updates the waitbar figure (if more than one solution file)
    if nExp > 1
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
        if h.Update(1,wStrNw,i/(1+nExp))
            [plotD,ok] = deal([],false);
            return
        else
            h.Update(2,'Calculating Scaled Coordinates',0);
        end
    end
    
    % sets the current experiment apparatus indices
    iApp = find(~cellfun('isempty',snTot(i).iMov.flyok));
    
    % sets the relevant x/y-locations for the current experiment    
    [dPx,dPy,R] = get2DCoordsBG(snTot(i),iApp);
    if (length(iApp) == 1)
        [dPx,dPy,R] = deal({dPx},{dPy},{R}); 
    end
                    
    % calculates the region exploration ratios for each grouping
    nApp2 = length(iApp);
    for k = 1:nApp2
        % updates the waitbar figure
        wStrNw = sprintf(...
                      '%s (Group %i of %i)',wStr{1+(nExp>1)},k,nApp2);
        if h.Update(1+(nExp>1),wStrNw,k/nApp2)
            [plotD,ok] = deal([],false);
            return
        end        
                
        % sets the new integer x/y coordinates (based on setup shape)
        switch mShape
            case 'Circle'
                % case is circular arenas
                RR = R{k}(:)';
                Px = num2cell(roundP(dPx{k}+RR),1);
                Py = num2cell(roundP(dPy{k}+RR),1);
                
            case 'Rectangle'
                % case is rectangular arenas                
                Px = num2cell(roundP(dPx{k}+R{k}(1,:)/2),1);
                Py = num2cell(roundP(dPy{k}+R{k}(2,:)/2),1);
        end
        
        % sets/appends the position look up table 
        DL = setupLookupTable(Px,Py,DL);
        
        % calculates the new exploration values
        switch mShape
            case 'Circle'        
                % case is circular arenas                
                Enw = cellfun(@(x,y,z)(detFlyExploration...
                    (x,y,DL,z,mShape)),Px,Py,num2cell(R{k}(:)'));
                
            case 'Rectangle'
                % case is rectangular arenas
                Enw = cellfun(@(x,y,z)(detFlyExploration...
                    (x,y,DL,z,mShape)),Px,Py,num2cell(R{k},1));
        end

%         ifok = snTot(i).iMov.flyok{j};
        ifok = 1:length(Enw);
        plotD(iApp(k)).E(1,ifok,i) = num2cell(Enw);       
    end
end

% sets trace values for each apparatus
for i = 1:nApp
    plotD(i) = calcMetricStats(plotD(i),'E');
end

% closes the waitbar figure
if ~h.Update(1,'Exploration Calculations Complete!',1)
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
[p,pP.pType,pP.plotErr] = deal(plotD{1}(ind),'BoxPlot',true);

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% retrieves the formatting struct
pF = retFormatStruct(pF,1);
xi = 1:length(ind);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% memory allocation
hAx = createSubPlotAxes(hP); 
hold(hAx,'on')

% sets the plot values for the 
plotBarBoxMetrics(hAx,xi,p,'E',pP,[],'b',100); 
        
% removes all the text labels
delete(findall(hAx,'type','text'))

% sets the plot
xTickLbl = snTot(1).iMov.pInfo.gName(ind)';
set(hAx,'ylim',[0 100],'xtick',xi,'xticklabel',xTickLbl)
  
% adds in the gridlines (if checked)
if (pP.plotGrid); grid(hAx,'on'); end

% ------------------------------ %
% --- PLOT AXES REFORMATTING --- %
% ------------------------------ %

% sets the x-axis labels
formatPlotAxis(hAx,pF,1);

% formats and resets the axis positions
resetAxesPos(hAx,1,1); 
setGroupString(hAx,pF,xi,snTot(1).iMov.pInfo.gName(ind),30);

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% ----------------------------------------------------------------------- %
% ---                         OTHER FUNCTIONS                         --- %
% ----------------------------------------------------------------------- %

%--- calculates the exploration ratio of the fly over the 2D arena
function E = detFlyExploration(x,y,DL,R,mShape)

% converts the x/y coordinates to cell arrays
[dR,X,Y] = deal(2,num2cell(x),num2cell(y));

switch mShape
    case 'Circle'
        % case is circular regions
        R = ceil(R);
        
    case 'Rectangle'
        % case is rectangular regions
        R = ceil(R) - 1;             
end

% calculates the inter-frame distance, and from this determine the frames
% where the distance is >= 2 pixels
[dX,dY] = deal(diff(x),diff(y));
D = sqrt(dX.^2 + dY.^2);
ii = find(D >= 2);

% for all frames where the interframe distance >1 pixel, set the paths
% between the points
if (~isempty(ii))
    % removes the small distance movements from the difference array
    [dX,dY] = deal(dX(ii),dY(ii));
    [sdX,sdY] = deal(num2cell(sign(dX)),num2cell(sign(dY)));
    
    % determines the frames where the differences are within the tables
    % limits, and resets the coords from here
    kk = (abs(dX) < size(DL,2)) & (abs(dY) < size(DL,1));
    ind = sub2ind(size(DL),abs(dY(kk))+1,abs(dX(kk))+1);
    
    % resets the coordinates
    X(ii(kk)) = cellfun(@(x,y,z)(x+z*y(:,1)),X(ii(kk)),DL(ind),sdX(kk),'un',0);
    Y(ii(kk)) = cellfun(@(x,y,z)(x+z*y(:,2)),Y(ii(kk)),DL(ind),sdY(kk),'un',0);
    
    % sets the very large interframe travelling flies and interpolates
    % their coordinates
    Xi = num2cell([x(ii(~kk)),x(ii(~kk)+1)],2);
    Yi = num2cell([y(ii(~kk)),y(ii(~kk)+1)],2);            
    
    % resets the x/y locations
    P = cellfun(@(x,y)(interpInterFrameCoords(x,y)),Xi,Yi,'un',0);
    X(ii(~kk)) = cellfun(@(x)(x(:,1)),P,'un',0);
    Y(ii(~kk)) = cellfun(@(x)(x(:,2)),P,'un',0);    
end

% ensures the indices are greater than zero
[X,Y] = deal(num2cell(max(cell2mat(X),1)),num2cell(max(cell2mat(Y),1)));

% calculates the number of pixels visited divided by the area of the
% circle (this gives the exploration ratio)
switch mShape
    case 'Circle'
        % case is circular regions
        XX = min(max(1,cell2mat(X)),2*R);
        YY = min(max(1,cell2mat(Y)),2*R);
        Apix = length(unique(sub2ind(2*R*[1 1],YY,XX)));
        E = min(1,Apix/(pi*(R-dR)^2));

    case 'Rectangle'
        % case is rectangular regions
        XX = min(max(1,cell2mat(X)),R(1));
        YY = min(max(1,cell2mat(Y)),R(2));
        Apix = length(unique(sub2ind(flip(R),YY,XX)));
        E = min(1,Apix/prod(R));
end
    
%
function DL = setupLookupTable(Px,Py,DLold)

% determines the most likely x/y coordinate dimension differences
[szX,szY] = detLikelyDim(Px,Py);

% memory allocation
DL = cell(szY+1,szX+1);
if (~isempty(DLold))
    % sets the old array within the current
    DL(1:size(DLold,1),1:size(DLold,2)) = DLold;
end

% interpolates the coordinates for each dimension difference
for i = 0:max(szY,(size(DLold,1)-1))
    for j = 0:max(szX,(size(DLold,2)-1))
        if (isempty(DL{i+1,j+1}))
            DL{i+1,j+1} = interpInterFrameCoords([0 j],[0 i]);
        end
    end
end

% --- determines the x/y-location differences that encapsulate the
%     proportion pTol of the entire location differences
function [szX,szY] = detLikelyDim(Px,Py)

% determines the 
pTol = 0.9999;

% calculates the difference in the x/y locations
dX = cell2mat(cellfun(@(x)(abs(diff(x))),Px,'un',0)');
dY = cell2mat(cellfun(@(x)(abs(diff(x))),Py,'un',0)');

% determines the x/y differences that encapsulate pTol proportion of the
% entire population of coordinate differences
[szX,szY] = deal(ceil(getThreshTol(dX,pTol)),ceil(getThreshTol(dY,pTol)));

% --- interpolates the x/y-coordinates between frames
function Pi = interpInterFrameCoords(X,Y)
    
% sets up the interpolation vector
Di = sqrt(diff(X).^2 + diff(Y).^2);
if (Di == 0); Pi = [X,Y]; return; end

% calculates the unique interpolated x/y pixel values
xi = linspace(0,1,ceil(1.5*Di))'; 
Pnw = roundP([interp1([0 1],X,xi),interp1([0 1],Y,xi)]);

% determines the unique values and removes the last line
[~,iA,~] = unique(Pnw,'rows');
[~,ii] = sort(iA);
Pi = Pnw(iA(ii(1:(end-1))),:);
