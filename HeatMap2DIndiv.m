% --- plot the individual 2D heat maps for a given experiment (2D analysis 
%     only)
function pData = HeatMap2DIndiv(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);

% sets the function name/
pData.Name = 'Individual 2D Heatmaps';
pData.Type = {'Indiv'};
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
    
    % sets the apparatus name/count
    pData.appName = snTotL.iMov.pInfo.gName;
    pData.nApp = length(pData.appName);
    
    % special parameter fields/data struct
%     [pData.hasSP,pData.hasRC,pData.hasRS] = deal(true,true,false);
    [pData.hasSR] = deal(true);
    pData.sP = initSpecialPara(snTotL,pData,[],1);
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

% determines if long experiment
isLong = ~detIfShortExpt(field2cell(snTot,'T'));

% initialises the parameter struct
nPara = 1 + isLong;
cP = setParaFields(nPara);

% sets the calculation parameter fields into the data struct
cP(1) = setParaFields([],'Number',25,'nGrid','Grid Resolution',[2 inf true]);
if (isLong)
    cP(2) = setParaFields([],'Boolean',false,'isDN','Separate Day/Night Activity');
end

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 3;
pP = setParaFields(nPara);

% sets the plot parameter fields into the data struct
pP(1) = setParaFields([],'Boolean',1,'pltLog','Plot Logarithmic Heatmap Values');
pP(2) = setParaFields([],'Number',5,'pDel','Heatmap Separation Gap',[1 20 true]);
pP(3) = setParaFields([],'Number',10,'nPlot','Max Sub-Region Heatmaps',[1 100 true]);

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
nApp = length(snTot.iMov.ok);
pF = setFormatFields(nApp);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],[],1);
pF.yLabel = setFormatFields([],[],1);
pF.Axis = setFormatFields([],[]);

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Title(i).String = snTot.iMov.pInfo.gName{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot,false,false);

% sets the output data parameter fields
oP = addYVarField(oP,'Heatmap','Ihm',[],7);

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
snTotL = snTot;
ok = checkFuncPara({'HasSubRegionStruct'},cP,snTotL);
if (~ok); plotD = []; return; end

% sets the day/night separation flag
if (isfield(cP,'isDN'))
    % flag is included in calculation parameters
    isDN = cP.isDN;
else
    % flag is not included in calculation parameters
    isDN = false;
end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensions
[ok,dnDel] = deal(true,5);

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,'Ihm',[]);

% other memory allocations
nG = cP.nGrid;
mShape = snTot.iMov.pInfo.mShape;

% determines the points outside of the circle
[y,x] = meshgrid(((1:nG)-0.5)/nG - 0.5); 
isN = sqrt(x.^2 + y.^2) > 0.5;
if (isDN); isN = [isN,true(nG,dnDel),isN]; end

% calculates the video frame rate and experiment apparatus indices
iApp = find(~cellfun('isempty',snTot.iMov.flyok));
szHM = nG*[1 (1+isDN)] + [0 dnDel*isDN];

% --------------------------------- %
% --- POSITION BIN CALCULATIONS --- %
% --------------------------------- %

% creates the waitbar figure
wStr = {'Retrieving 2D Coordinates'};
h = ProgBar(wStr,'2D Heatmap Calculations');

% retrieves the 2D coordinates (based on region shape)
[dPx,dPy,R] = get2DCoordsBG(snTot,iApp);

% ensures the data is stored in cell arrays
if ~iscell(dPx)
    [dPx,dPy,R] = deal({dPx},{dPy},{R});
end

% determines the day/night time points
if isDN
    isDay = detDayNightTimePoints(snTot);
else
    isDay = true(size(dPx{1},1),1);
end

% determines which trace belongs to which sub-region. from this, scale
% the data values
for j = 1:length(iApp)
    % updates the waitbar figure (if more than one solution file)
    wStrNw = sprintf('Calculating Heatmaps (Region %i of %i)',...
                        j,length(iApp));
    if h.Update(1,wStrNw,0.5*(1+j/(1+length(iApp))))
        [plotD,ok] = deal([],false);
        return
    end
    
    % memory allocation
    [Px,Py] = deal(cell(1,1+isDN));
    
    % sets the x/y direction scale factors (based on shape)
    switch mShape
        case 'Circle'
            % case is circular regions
            [xScl,yScl] = deal(2*R{j}(:)');
            
        case 'Rectangle'
            % case is rectangular regions
            [xScl,yScl] = deal(R{j}(1,:),R{j}(2,:));            
    end
    
    % calculates the normalised x-locations
    for k = 1:(1+isDN)
        % sets the day/night flags for the current phase
        if (k == 1)
            iDN = isDay;
        else
            iDN = ~isDay;
        end    
    
        % calculates the normalised coordinates
        Px{k} = num2cell(dPx{j}(iDN,:)./xScl + 1/2,1);
        Py{k} = num2cell(dPy{j}(iDN,:)./yScl + 1/2,1); 
    end
    
    % sets the heatmap
    Ihm = cell(size(dPx{j},2),1);
    for i = 1:length(Ihm)
        % memory allocation
        Ihm{i} = zeros(szHM);
        
        % sets the individual heatmaps
        for k = 1:(1+isDN)   
            % sets the column indices
            iCol = (k-1)*(nG+dnDel) + (1:nG); 
            
            % sets the x/y coordinates of the points
            iPx = max(1,min(roundP(Px{k}{i}*(nG-1))+1,nG));
            iPy = max(1,min(roundP(Py{k}{i}*(nG-1))+1,nG));   

            % sets the heatmap array for the given group
            if (~isempty(iPx))
                ind = sub2ind(nG*[1 1],iPy,iPx);
                Ihm{i}(:,iCol) = Ihm{i}(:,iCol) + ...
                                    reshape(hist(ind,1:(nG^2)),nG*[1 1]);            
            end
        end
        
        % removes the outside region values
        Ihm{i}(isN) = NaN;        
    end
    
    % sets the heatmap array into the plot data struct
    plotD(iApp(j)).Ihm = Ihm;
end

% closes the waitbar figure
if ~h.Update(1,'2D Heatmap Calculations Complete!',1)
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

% memory allocation and other dimensioning
ind = sP.pInd;
p = plotD{1}(ind);
mShape = snTot.iMov.pInfo.mShape;

% if there are no heatmaps for this group, then exit
if (isempty(p.Ihm)); return; end

% retrieves the heatmap array
Ihm = p.Ihm;
if (pP.pltLog)
    Ihm = cellfun(@(y)(log10(y+1)),Ihm,'un',0);
end

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');
pPos = get(hP,'Position');

% reduces the output plots to the required amount
nHm = length(Ihm);
Ihm = Ihm(sort(randperm(nHm,min(nHm,pP.nPlot))));

% sets the custom colormap values
cMap = [zeros(1,3);colormap('jet');ones(1,3)];

% determines the maximum values over all the heatmaps
Ymx = max(cellfun(@(x)(max(x(:))),Ihm));
dcMap = 1/(size(cMap,1)-3);

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% resets the title strings
pF.Title = pF.Title(ind);

% retrieves the formatting struct
pF = retFormatStruct(pF,1);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% determines the dimensions of the plot grid
nR = floor(sqrt(pP.nPlot*pPos(4)/pPos(3)));
nC = ceil(pP.nPlot/nR);
[szR,szC] = size(Ihm{1});

% creates the subplot axis
hAx = createSubPlotAxes(hP,[1,1],1);  
axis(hAx,'off')

% memory allocations   
IhmT = NaN(nR*szR+(nR+1)*pP.pDel,nC*szC+(nC+1)*pP.pDel);
   
% sets the individual heatmaps into the combined heatmap
for j = 1:length(Ihm)
    % sets the global row/column indices
    [iRow,iCol] = deal(floor((j-1)/nC)+1,mod(j-1,nC)+1);
    
    % sets the row/column indices for the new heatmap
    iR = iRow*pP.pDel + (iRow-1)*szR + (szR:-1:1);
    iC = iCol*pP.pDel + (iCol-1)*szC + (1:szC);    

    % sets the individual heatmap
    IhmT(iR,iC) = Ihm{j};
end

% reverses the order of the arrays
IhmT = IhmT(end:-1:1,:);

% sets the heatmap edge and outline colours
switch mShape
    case 'Circle'
        % case is circular regions
        Bnan = isnan(IhmT);
        IhmT(bwmorph(~Bnan,'dilate',1) & Bnan) = -dcMap;
        IhmT(isnan(IhmT)) = Ymx + dcMap;
        
    case 'Rectangle'
        % case is rectangular regions
        IhmT(isnan(IhmT)) = 0;
end

% creates the heatmap image    
set(get(hP,'parent'),'CurrentAxes',hAx);
imagesc(IhmT);
colormap(cMap)    

% ensures the axis is equal is size and hold is on
axis(hAx,'equal');         
set(hAx,'xtick',[],'ytick',[],'yticklabel',[],'xticklabel',[])     
drawnow; 

% formats the plot axis
formatPlotAxis(hAx,pF,1);      

% sets the x/y axis limits    
set(hAx,'box','off','xcolor','w','ycolor','w')
xL = get(hAx,'xlim');
set(hAx,'xlim',xL + 0.001*diff(xL)*[-1 1]);

% turns the axis on
axis(hAx,'on')

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)
