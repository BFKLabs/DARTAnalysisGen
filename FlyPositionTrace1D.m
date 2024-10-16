% --- plots the x-location of the fly over the duration of the experiment 
%     (1D analysis only)
function pData = FlyPositionTrace1D(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = '1D Fly Position Traces';
pData.Type = {'Indiv'};
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
    [pData.hasTime,pData.hasSR] = deal(true);    
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
rI.Dur = 'None';
rI.Shape = '1D';
rI.Stim = 'None';
rI.Spec = 'None';
        
% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 2;
cP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the plot parameter fields into the data struct
cP(1) = setParaFields(a{1},'Boolean',true,'useDS','Use Signal Downsampling');
cP(2) = setParaFields(a{1},'Number',50,'sRate','Signal Downsample Rate',[10 1000 1],{1,2});

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 3;                       
pP = setParaFields(nPara);
isLong = ~detIfShortExpt(field2cell(snTot,'T'));

% sets the tab list names
a = {'1 - General'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},'Boolean',isLong,'pltDN','Plot Day/Night Background Image',[],{0,~isLong});
pP(2) = setParaFields(a{1},'Boolean',0,'isZeitG','Use Zeitgeiber Time Format');
pP(3) = setParaFields(a{1},'Boolean',1,'pltRej','Plot Rejected Flies');

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
nApp = length(snTot.iMov.ok);  
pF = setFormatFields(1);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],'Time');
pF.yLabel = setFormatFields([],'Fly Index');
pF.Axis = setFormatFields([],[]);

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Title(i).String = snTot.iMov.pInfo.gName{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot,false,false,false);

% sets the independent variable fields
oP = addXVarField(oP,'Time','T','Time');

% sets the dependent variable fields
oP = addYVarField(oP,'Position','YD',[],5,{'T'},1);

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

% other initialisations
flyok = snTot.iMov.flyok;
[nApp,ok] = deal(length(flyok),true);
T = cell2mat(snTot.T);

% sets the down-sampling rate indices
if cP.useDS
    ii = 1:cP.sRate:length(T);
else
    ii = 1:length(T);
end

% memory allocation
plotD = initPlotValueStruct(snTot,pData,cP,...
                                'T',T(ii),'Y',[],...
                                'TD',[],'YD',[]);

% memory allocation                            
nDay = size(plotD(1).YD,1);  
szR = [snTot.iMov.nCol,snTot.iMov.nRow];

% ---------------------------- %
% --- FLY LOCATION SETTING --- %
% ---------------------------- %

% determines the 
iC = snTot.iMov.iC;
sFac = snTot.sgP.sFac;

% sets the position values for all the apparatus
for i = 1:nApp
    % 
    indR = sub2ind(szR,snTot.cID{i}(:,2),snTot.cID{i}(:,1));
    xOfs = sFac*arrayfun(@(x)(iC{x}(1)-1),indR)'; 
    xRng = sFac*arrayfun(@(x)(range(iC{x})),indR)';
    Px = (snTot.Px{i}-xOfs)./xRng;
        
    % retrieves the array of x-locations
    plotD(i).Y = 1 - Px(ii,:);    
    xPx = snTot.Px{i}(ii,:);
    
    % sets the speed values for each of the days    
    plotD(i).YD = repmat({NaN(length(ii),1)},1,size(plotD(i).YD,2));
    for j = 1:nDay
        for k = 1:size(xPx,2)
            plotD(i).YD{j,k,1} = xPx(:,k);
        end
    end
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

% memory allocation and array indexing
ind = sP.pInd;
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% sets the plot indices
ind0 = find(p.T >= sP.xLim(1),1,'first');
indF = find(p.T <= sP.xLim(2),1,'last');

% resets the plot arrays
[tMlt,tUnits] = getTimeScale(snTot.T{end}(end));

% sets the other parameters
flyok = snTot.iMov.flyok{ind};
[nFly,Px] = deal(size(p.Y,2),p.Y);

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the formatting data struct
pF = retFormatStruct(pF,1);

% resets the axis font size if greater than threshold
pF.Axis(1).Font.FontSize = pF.Axis(1).Font.FontSize - 4*(nFly > 20);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% plots the day/night bands (if required)
if pP.pltDN
    hAx = plotDayNightGraph(hP,snTot,p.T*tMlt,nFly+1,1,[1,1],tMlt);
else
    hAx = createSubPlotAxes(hP);
    set(hAx,'xlim',p.T([1 end]));
end

% if the plots do not exist, then create them for each apparatus
hold(hAx,'on'); axis(hAx,'on')
for i = 1:nFly
    % sets the plot colour/position values
    if flyok(i)
        pCol = getTraceColour(mod(i-1,4)+1);        
    else
        if pP.pltRej
            pCol = 'k';
        else
            pCol = [];
        end
    end
    
    % calculates the normalised positions
    if ~isempty(pCol)
        PxNw = Px(:,i) + ((i-1) + 0.5);

        % plots the trace
        plot(hAx,p.T*tMlt,PxNw,'color',pCol,'LineWidth',0.5);                
        if (i ~= nFly)
            % adds the seperator markers
            plot(get(hAx,'xlim'),(i+0.5)*[1 1],'k','linewidth',1)
        end
    else
        [xL,yL] = deal(get(hAx,'xlim'),(i-1)+[0 1]+0.5);
        fill(xL([1 2 2 1]),yL([1 1 2 2]),'k','facealpha',0.3);
    end
end    

%
% Tadd = floor(convertTime(T(1),'hrs','sec'));
% snTot.iExpt.Timing.T0 = datevec(addtodate(datenum(snTot.iExpt.Timing.T0),Tadd,'s'));
    
% case is using the absolute time axis
if pP.isZeitG            
    % case is using Zeitgeiber time
    setZeitGTimeAxis(hAx,p.T,snTot);   
    pF.xLabel.String = 'Zeitgeiber Time';
else
    % case is not using Zeitgeiber time    
    if ((pP.pltDN) || (strcmp(tUnits,'Hours')))
        % if plotting day/night, then set absolute time         
        setAbsTimeAxis(hAx,p.T,snTot);    
        pF.xLabel.String = 'Time';
    else        
        pF.xLabel.String = 'Time (min)';
    end
end

% sets the axis properties
axis(hAx,'ij')
xLim = p.T([ind0,indF])*tMlt;
set(hAx,'xLim',xLim,'ylim',[1 nFly]+0.5*[-1 1],'ytick',1:nFly,...
        'yticklabel',num2cell(1:nFly),'box','on','UserData',1)

% resets the axis positions
formatPlotAxis(hAx,pF,sP.pInd);
resetAxesPos(hAx,1,1,[0.025,0]); 

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)
