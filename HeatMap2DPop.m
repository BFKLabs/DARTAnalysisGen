% --- plot the combined population 2D heat maps (2D analysis only)
function pData = HeatMap2DPop(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);

% sets the function name/
pData.Name = 'Population 2D Heatmaps';
pData.Type = {'Pop','Multi'};
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
    [pData.hasSP,pData.hasRC,pData.hasRS] = deal(true,true,false);
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

% determines if long experiment
isLong = ~detIfShortExpt(field2cell(snTot,'T'));

% initialises the parameter struct
nPara = 1 + isLong;
cP = setParaFields(nPara);

% sets the calculation parameter fields into the data struct
cP(1) = setParaFields([],'Number',25,'nGrid','Grid Resolution',[2 inf true]);
if isLong
    cP(2) = setParaFields([],'Boolean',false,'isDN','Separate Day/Night Activity');
end

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 1;
pP = setParaFields(nPara);

% sets the plot parameter fields into the data struct
pP(1) = setParaFields([],'Boolean',1,'pltLog','Plot Logarithmic Heatmap Values');

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
oP = addYVarField(oP,'Population Heatmap','Ihm',[],6);

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
if ~ok; plotD = []; return; end

% sets the day/night separation flag
if isfield(cP,'isDN')
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
dnDel = 5;
mShape = snTotL.iMov.pInfo.mShape;
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,'Ihm',[]);

% other memory allocations
nG = cP.nGrid;
[Px,Py] = deal(cell(1,nApp,1+isDN));

% determines the points outside of the circle
switch mShape
    case 'Circle'
        [y,x] = meshgrid(((1:nG)-0.5)/nG - 0.5); 
        isN = sqrt(x.^2 + y.^2) > 0.5;
        if isDN; isN = [isN,true(nG,dnDel),isN]; end
end

% memory allocation
szHM = nG*[1 (1+isDN)] + [0 dnDel*isDN];
Ihm = repmat({zeros(szHM)},1,nApp);

% --------------------------------- %
% --- POSITION BIN CALCULATIONS --- %
% --------------------------------- %

% creates the waitbar figure
wStr = {'Setting 2D Coordinates'};
h = ProgBar(wStr,'2D Heatmap Calculations');

% loops through each of the experiments calculating the velocity values
for i = 1:nExp
    % updates the waitbar figure (if more than one solution file)
    wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
    if h.Update(1,wStrNw,i/(1+nExp))
        [plotD,ok] = deal([],false);
        return
    end

    % calculates the video frame rate and experiment apparatus indices
    iApp = find(~cellfun('isempty',snTot(i).iMov.flyok));    
    
    % retrieves the 2D coordinates (based on region shape)
    [dPx,dPy,R] = get2DCoordsBG(snTot(i),iApp);    
        
    % ensures the data is stored in cell arrays
    if (length(iApp) == 1)
        [dPx,dPy,R] = deal({dPx},{dPy},{R});
    end
    
    % determines the day/night time points
    if isDN
        isDay = detDayNightTimePoints(snTot(i));
    else        
        isDay = true(size(dPx{1},1),1);  
    end        
    
    % determines which trace belongs to which sub-region. from this, scale
    % the data values
    nApp2 = length(iApp);
    for j = 1:nApp2
        % sets the x/y direction scale factors (based on shape)
        switch mShape
            case 'Circle'
                % case is circular regions
                [xScl,yScl] = deal(2*R{j}(:)');

            case 'Rectangle'
                % case is rectangular regions
                [xScl,yScl] = deal(R{j}(1,:),R{j}(2,:));            
        end        
        
        % sets the radius values       
        jj = iApp(j);        
        for k = 1:(1+isDN)
            % sets the day/night flags for the current phase
            if (k == 1)
                iDN = isDay;
            else
                iDN = ~isDay;
            end
            
            % sets the normalised x/y coordinates
            Px{1,jj,k} = [Px{1,jj,k},...
                    num2cell(dPx{j}(iDN,:)./xScl + 1/2,1)];
            Py{1,jj,k} = [Py{1,jj,k},...
                    num2cell(dPy{j}(iDN,:)./yScl + 1/2,1)];
        end
    end
end    

% sets the final combined heatmaps
for i = 1:nApp
    for k = 1:(1+isDN)
        iCol = (k-1)*(nG+dnDel) + (1:nG);        
        for j = 1:length(Px{1,i,k})
            % determines the x/y coordinate indices
            iPx = max(1,min(roundP(Px{1,i,k}{j}*(nG-1))+1,nG));
            iPy = max(1,min(roundP(Py{1,i,k}{j}*(nG-1))+1,nG));     

            % sets the heatmap array for the given group
            if ~isempty(iPx)
                ind = sub2ind(nG*[1 1],iPy,iPx);      
                dIhm = reshape(hist(ind,1:(nG^2)),nG*[1 1]);
                Ihm{i}(:,iCol) = Ihm{i}(:,iCol) + dIhm;
            end
        end
    end
    
    % removes the outside regions
    if exist('isN','var')    
        Ihm{i}(isN) = NaN;
    end
    
    % sets the heatmap array into the plot data struct
    plotD(i).Ihm = Ihm{i};   
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

% sets the plotting indices and subplot indices
[ind,m,n] = deal(find(sP.Sub.isPlot),sP.Sub.nRow,sP.Sub.nCol);
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% retrieves the heatmap array
Ihm = field2cell(plotD{1}(ind),'Ihm')';
if pP.pltLog
    Ihm = cellfun(@(x)(log10(x+1)),Ihm,'un',0);
end

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% memory allocation
hAx = cell(nApp,1);

% sets the custom colormap values
cMap = [zeros(1,3);colormap('jet');ones(1,3)];

% determines the maximum values over all the heatmaps
Ymx = max(cellfun(@(x)(max(x(:))),Ihm));
dcMap = 1/(size(cMap,1)-3);

% parameters
mShape = snTot(1).iMov.pInfo.mShape;

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% resets the title strings
pF.Title = pF.Title(ind);

% retrieves the formatting struct
pF = retFormatStruct(pF,nApp);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% loops through all the subplots
for i = 1:nApp
    % creates the subplot axis
    hAx{i} = createSubPlotAxes(hP,[m,n],i);
    
    % sets the heat-map into a larger array
    IhmF = NaN(size(Ihm{i}) + 2);
    IhmF(1+(1:size(Ihm{i},1)),1+(1:size(Ihm{i},2))) = Ihm{i};
    
    % sets the heatmap edge and outline colours
    switch mShape
        case 'Circle'
            % case is circular regions
            Bnan = isnan(IhmF);        
            IhmF(~Bnan & Bnan) = -dcMap;
            IhmF(isnan(IhmF)) = Ymx + dcMap;
            
        case 'Rectangle'
            % case is rectangular regions
            IhmF(isnan(IhmF)) = Ymx + dcMap;                        
    end
    
    % creates the heatmap image
    imagesc(IhmF);
    colormap(cMap);   
    set(hAx{i},'clim',[-dcMap,(Ymx+dcMap)])
    
    % ensures the axis is equal is size and hold is on
    axis(hAx{i},'equal');         
    set(hAx{i},'xtick',[],'ytick',[],'yticklabel',[],'xticklabel',[])     
    
    % formats the plot axis
    formatPlotAxis(hAx{i},pF,i);      
        
    % sets the x/y axis limits    
    set(hAx{i},'box','off','xcolor','w','ycolor','w')
    xL = get(hAx{i},'xlim');
    set(hAx{i},'xlim',xL + 0.001*diff(xL)*[-1 1]);
end    
    
% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)
