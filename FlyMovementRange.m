% --- displays the aggregate distance travelled each minute by a population
%     of flies over the duration of an experiment(s)
function pData = FlyMovementRange(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Fly Distance Statistics';
pData.Type = {'Pop','Multi'};
pData.fType = [3 1 1 1];
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
    [pData.hasSR,pData.hasRC] = deal(true,false);    
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
rI.Shape = 'None';
rI.Stim = 'None';
rI.Spec = 'None';
        
% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 1;
cP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the plot parameter fields into the data struct
cP(1) = setParaFields(a{1},'Number',1,'Bsz','Histogram Bin Size (mm)',[1 10 true]);

% sets the tool-tip strings
cP(1).TTstr = 'Width of the bin size used to calculated the histogram';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 1;                       
pP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},'Number',0.95,'pCDF','CDF Percentage Marker',[0.01 0.99 false]);

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
nPlot = 2;
pF = setFormatFields(nPlot);

% sets the font structs
tpF = setupFontStruct('FontSize',18);
lpF = setupFontStruct('FontSize',16);
apF = setupFontStruct('FontSize',14);

% initialises the font structs
pF.Title(1) = setFormatFields(tpF,'Distance Histogram');
pF.xLabel(1) = setFormatFields(lpF,'Distance (mm/min)');
pF.yLabel(1) = setFormatFields(lpF,'Frequency');
pF.Axis(1) = setFormatFields(apF,[]);

% initialises the font structs
pF.Title(2) = setFormatFields(tpF,'Distance Empirical CDF');
pF.xLabel(2) = setFormatFields(lpF,'Distance (mm/min)');
pF.yLabel(2) = setFormatFields(lpF,'Cumulative %age');
pF.Axis(2) = setFormatFields(apF,[]);

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot,false,false,false,true);

% sets the independent output variables
oP = addXVarField(oP,'Freq Dist (mm)','XH','Other');
oP = addXVarField(oP,'Prop Dist (mm)','xCDF','Other');

% sets the dependent output variables
oP = addYVarField(oP,'Frequency','NH',[],4,{'XH'});
oP = addYVarField(oP,'Proportion','fCDF',[],4,{'xCDF'});

% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% field retrievals
iAx = dcObj.getSelectAxesIndex;

% sets the common class fields
dcObj.xName = 'Distance';
dcObj.xUnits = 'mm/min';
dcObj.grpName = dcObj.pData.appName;

% sets up the metric specific class fields
switch iAx
    case 1
        % case is the distance histogram
        dcObj.yName = 'Frequency';
        dcObj.yUnits = 'Count';
        dcObj.pType = 'Bar Graph';
        dcObj.combFig =  false;
        dcObj.xGrp = dcObj.plotD{1}(1).XH;
        
    case 2
        % case is the distance cdf
        dcObj.yName = 'Cumulative %';
        dcObj.yUnits = '%';
        dcObj.pType = 'Trace';
        dcObj.combFig =  true;
        
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
    cP = retParaStruct(pData.cP,gPara);  
end

% sets the movement type to absolute range
cP.movType = 'Absolute Distance';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
[nApp,ok] = deal(length(snTot(1).iMov.ok),true);
nExp = length(snTot);

% sets the calculation parameters
DMax = 1000;
[cP.tBin,cP.nGrp] = deal(60,'1');

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                                'NH',[],'XH',[],'xCDF',[],'fCDF',[]);

% --------------------------- %
% --- METRIC CALCULATIONS --- %
% --------------------------- %

% memory allocation
V = cell(nApp,nExp);

% creates the waitbar figure
wStr = {'Overall Progress','Movement CDF Calculations'};
wOfs = (nExp > 1); wStr = wStr((2-wOfs):end);
h = ProgBar(wStr,'Movement Range Calculations');

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

    % sets the total time and ok flags
    Ttot = cell2mat(snTot(i).T);
    flyok = snTot(i).iMov.flyok;

    % determines the binned indices (for time length, tBin) and determines 
    % the bins which has at least two time points
    h.Update(1+wOfs,'Determining Time Bins',0.50);
    indB = detTimeBinIndices(Ttot,cP.tBin);
    indB(cellfun('length',indB) < cP.tBin/2) = {[]};
    
    % calculates the sleep metrics
    for j = 1:nApp       
        % updates the waitbar figure
        wStrNw = sprintf('Calculating Fly Range (Region %i of %i)',j,nApp);
        if h.Update(1+wOfs,wStrNw,0.50*(1+(j/(nApp+1))))
            [plotD,ok] = deal([],false);
            return            
        end
         
        % only calculate if values exist...
        if ~isempty(snTot(i).Px{j})       
            % sets the wake metrics
            [~,~,~,V{j,i}] = ...
                    calcWakingMetrics(snTot(i),Ttot,indB,cP,j,flyok{j}); 
        end
    end
    
    % updates the waitbar figure
    h.Update(1+wOfs,'Fly Range Calculations Complete!',1); 
end

% sets the values into the plot data struct
for i = 1:nApp    
    % determines the non-empty arrays
    ii = ~cellfun('isempty',V(i,:));
    
    % converts the total velocity array to a vector
    VT = combineNumericCells(cellfun(@(x)(cell2mat(x)),num2cell(...
        cellfun(@(x)(cell2mat(x')'),V(i,ii),'un',0),1),'un',0));
    VT = VT(VT<DMax);
    
    % calculates the histogram/ecdf values
    [NH,XH] = hist(VT,(0.00:cP.Bsz:(DMax-cP.Bsz))+cP.Bsz/2); 
    [plotD(i).NH,plotD(i).XH] = deal(NH(1:end)',XH(1:end)');
    
    % calculates the ecdf values (then interpolate to fixed values)
    if isempty(VT)
        [xCDF,fCDF] = deal([0,DMax],[0,0]);
    else
        [fCDF,xCDF] = ecdf(VT);         
    end
    
    if length(xCDF) == 2
        [plotD(i).xCDF,xCDF] = deal([0,DMax]);
        plotD(i).fCDF = interp1(xCDF,[0,0],xCDF);
    else
        plotD(i).xCDF = linspace(xCDF(1),xCDF(end),1001)';
        plotD(i).fCDF = interp1(xCDF(2:end),fCDF(2:end),plotD(i).xCDF);
    end
end

% closes the waitbar figure
h.closeProgBar(); 

% ----------------------------------------------------------------------- %
% ---                        PLOTTING FUNCTION                        --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,varargout] = plotFunc(snTot,pData,plotD,ind)

% retrieves the plotting paraeter struct
pP = retParaStruct(pData.pP);
cP = retParaStruct(pData.cP);
sP = retParaStruct(pData.sP);
pF = pData.pF;

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the plot values data struct
p = plotD{1};

% array indexing
[ind,nApp] = deal(sP.pInd,length(p));
hPlt = cell(nApp,1);

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% resets the formatting struct
pF = retFormatStruct(pF,2);

% sets the legend strings
grpStr = snTot(1).iMov.pInfo.gName{ind};
pF.Legend.String = snTot(1).iMov.pInfo.gName;
pF.Title(1).String = sprintf('%s (%s)',pF.Title(1).String,grpStr);

% sets the colours
col = num2cell(distinguishable_colors(nApp,'w'),2);
yLim = [1 10^(ceil(log10(max(p(ind).NH(2:end)))))];
if (yLim(2) == 0) || (range(yLim) == 0)
    yLim(2) = 10;
end

% sets the x-axis limits
[xCDF,fCDF] = field2cell(p,{'xCDF','fCDF'});
if any(cellfun(@range,fCDF)>0)
    xLim = [0 ceil(max(cellfun(@(x,y)...
                (x(find(y>0.999,1,'first'))),xCDF,fCDF))/50)*50];
else
    xLim = [0 max(cellfun(@(x)(x(end)),xCDF))];
end

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');             
                
% memory allocation
hAx = cell(2,1);

% ---------------------------------- %
% --- DISTANCE HISTOGRAM SUBPLOT --- %
% ---------------------------------- %

% plots the empirical cdf
hAx{1} = createSubPlotAxes(hP,[2,1],1); hold on; 

% sets the histogram plot range
xMx = ceil(find(p(ind).NH>0,1,'last')/(100/cP.Bsz))*100;
iMx = 2:find((p(ind).XH <= xMx),1,'last');

% plots the histogram
bar(hAx{1},p(ind).XH(iMx),p(ind).NH(iMx)+1); 
set(hAx{1},'ylim',yLim,'yscale','log','ygrid','on',...
           'xgrid','on','UserData',1)

% sets the plot labels
formatPlotAxis(hAx{1},pF,1);    

% ---------------------------- %
% --- DISTANCE CDF SUBPLOT --- %
% ---------------------------- %

%
hAx{2} = createSubPlotAxes(hP,[2,1],2); hold on; 

% plots the CDF values for each of the apparatus
for i = 1:nApp
    % plots the CDF values
    stairs(hAx{2},p(i).xCDF,100*p(i).fCDF,'color',col{i},'UserData',i); 

    % plots the other markers
    xx = p(i).xCDF(find(p(i).fCDF <= pP.pCDF,1,'last'));
    hPlt{i} = plot(hAx{2},xx*[1 1],100*[0 pP.pCDF],'--',...
        'linewidth',2,'color',col{i},'UserData',i,'HitTest','off'); 
end
   
% creates the legend object
hLg = createLegendObj(hPlt,pF.Legend);

% sets the plot labels
set(hAx{2},'ylim',[0 100],'xlim',xLim,'box','on','xgrid','on',...
           'ygrid','on','UserData',2)
formatPlotAxis(hAx{2},pF,2);   

% formats and resets the axis positions
resetAxesPos(hAx,2,1,[0.0625 0]); 
resetLegendAxisPos(hAx,hLg,[2 1])

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

a = 1;
