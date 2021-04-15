% --- plot the individual 2D heat maps for a given experiment (2D analysis 
%     only)
function pData = HeatMap2DIndiv(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);

% sets the function name/
pData.Name = 'Individual 2D Heatmaps';
pData.Type = 'Indiv';
pData.fType = [2 1 1 3];

% initialises the other fields  (if input argument provided)
if (nargin == 1)
    % parameter data struct initialisation
    snTotL = snTot(1);
    pData.cP = initCalcPara(snTot);
    pData.pP = initPlotPara(snTot);
    pData.oP = initOutputPara(snTot);
    pData.pF = initPlotFormat(snTotL);
    
    % sets the apparatus name/count
    pData.appName = snTotL.appPara.Name;
    pData.nApp = length(snTotL.appPara.Name);
    
    % special parameter fields/data struct
    [pData.hasSP,pData.hasRC,pData.hasRS] = deal(true,true,false);
    pData.sP = initSpecialPara(snTotL,pData);
end

% ----------------------------------------------------------------------- %
% ---                 PARAMETER STRUCT SETUP FUNCTIONS                --- %
% ----------------------------------------------------------------------- %

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
pP(3) = setParaFields([],'Number',10,'nPlot','Max Sub-Region Heatmaps',[1 20 true]);

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
nApp = length(snTot.appPara.ok);
pF = setFormatFields(nApp);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],[],1);
pF.yLabel = setFormatFields([],[],1);
pF.Axis = setFormatFields([],[]);

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Title(i).String = snTot.appPara.Name{i};
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
[nApp,ok,dnDel] = deal(length(snTot.appPara.ok),true,5);

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,'Ihm',[]);

% other memory allocations
nG = cP.nGrid;

% determines the points outside of the circle
[y,x] = meshgrid(((1:nG)-0.5)/nG - 0.5); 
isN = sqrt(x.^2 + y.^2) > 0.5;
if (isDN); isN = [isN,true(nG,dnDel),isN]; end

% calculates the video frame rate and experiment apparatus indices
iApp = find(~cellfun(@isempty,snTot.appPara.flyok));
szHM = nG*[1 (1+isDN)] + [0 dnDel*isDN];

% --------------------------------- %
% --- POSITION BIN CALCULATIONS --- %
% --------------------------------- %

% creates the waitbar figure
wStr = {'Retrieving 2D Coordinates'};
h = ProgBar(wStr,'2D Heatmap Calculations');

% sets the relevant x/y-locations for the current experiment  
[dPx,dPy,R] = get2DCoordsBG(snTot,iApp);

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
    
    % calculates the normalised x-locations
    RR = repmat(R{j}(:)',size(dPx{j},1),1);
    for k = 1:(1+isDN)
        % sets the day/night flags for the current phase
        if (k == 1)
            iDN = isDay;
        else
            iDN = ~isDay;
        end    
    
        Px{k} = num2cell((dPx{j}(iDN,:)+RR(iDN,:))./(2*RR(iDN,:)),1);
        Py{k} = num2cell((dPy{j}(iDN,:)+RR(iDN,:))./(2*RR(iDN,:)),1); 
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

% sets the plotting indices and subplot indices
[ind,m,n] = deal(find(sP.Sub.isPlot),sP.Sub.nRow,sP.Sub.nCol);
nApp = length(ind); if (nApp == 0); return; end

% retrieves the heatmap array
Ihm = field2cell(plotD{1}(ind),'Ihm')';
if (pP.pltLog)
    Ihm = cellfun(@(x)(cellfun(@(y)(log10(y+1)),x,'un',0)),Ihm,'un',0);
end

% reduces the output plots to the required amount
indH = cellfun(@(x)(sort(randperm(length(x),min(length(x),pP.nPlot)))),Ihm,'un',0);
Ihm = cellfun(@(x,y)(x(y)),Ihm,indH,'un',0);

% retrieves the parent axis
hP = get(gca,'parent');

% memory allocation
hAx = cell(nApp,1);

% sets the custom colormap values
cMap = [zeros(1,3);colormap('jet');ones(1,3)];

% determines the maximum values over all the heatmaps
Ymx = max(cellfun(@(x)(max(cellfun(@(y)(max(y(:))),x(:)))),Ihm));
dcMap = 1/(size(cMap,1)-3);

% paramters
nG = cP.nGrid;

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

% memory allocation
IhmT = cell(1,length(Ihm));

% 
szHm = cell2mat(cellfun(@(x)(size(x)),Ihm(:),'un',0));
[nR,nC] = deal(max(szHm(:,1)),max(szHm(:,2)));
[szR,szC] = size(Ihm{1}{1});

% sets the combine heatmap arrays
for i = 1:length(IhmT)
    % creates the subplot axis
    hAx{i} = createSubPlotAxes(hP,[m,n],i);  
    axis(hAx{i},'off')
    
    % memory allocations
    [nRnw,nCnw] = size(Ihm{i});    
    IhmT{i} = NaN(nR*szR+(nR+1)*pP.pDel,nC*szC+(nC+1)*pP.pDel);
   
    % sets the individual heatmaps into the combined heatmap
    for j = 1:nRnw
        for k = 1:nCnw
            % sets the row/column indices for the new heatmap
            iR = j*pP.pDel + (j-1)*szR + (szR:-1:1);
            iC = k*pP.pDel + (k-1)*szC + (1:szC);
            
            % sets the individual heatmap
            IhmT{i}(iR,iC) = Ihm{i}{j,k}; 
        end
    end
    
    % reverses the order of the arrays
    IhmT{i} = IhmT{i}(end:-1:1,:);
    
    % sets the overall maximum values
    Bnan = isnan(IhmT{i});
    IhmT{i}(bwmorph(~Bnan,'dilate',1) & Bnan) = -dcMap;
    IhmT{i}(isnan(IhmT{i})) = Ymx+dcMap;
end

% loops through all the subplots
for i = 1:nApp    
    % creates the heatmap image    
    set(get(hP,'parent'),'CurrentAxes',hAx{i});
    imagesc(IhmT{i});
    colormap(cMap)    
        
    % ensures the axis is equal is size and hold is on
    axis(hAx{i},'equal');         
    set(hAx{i},'xtick',[],'ytick',[],'yticklabel',[],'xticklabel',[])     
    if (nApp == 1); drawnow; end
    
    % formats the plot axis
    formatPlotAxis(hAx{i},pF,i);      
        
    % sets the x/y axis limits    
    set(hAx{i},'box','off','xcolor','w','ycolor','w')           
    if (~isHG1)
        xL = get(hAx{i},'xlim');
        set(hAx{i},'xlim',xL + 0.001*diff(xL)*[-1 1]);
    end    
    
    % turns the axis on
    axis(hAx{i},'on')
end

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)
