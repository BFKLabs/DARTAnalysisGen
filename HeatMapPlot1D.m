% --- plots the heatmap of the fly location over the duration of an 
%     experiment (1D experiment Only)
function pData = HeatMapPlot1D(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = '1D Fly Location Heatmap (Total)';
pData.Type = {'Pop'};
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
    [pData.hasSP,pData.hasRC,pData.hasRS] = deal(true,false,false);
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
rI.Shape = '1D';
rI.Stim = 'None';
rI.Spec = 'None';
        
% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 3;
cP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the parameter list for 
pList = {'Mean Location','Median Location'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'Number',60,'tBin','Time Bin Size (sec)',[5 3600 false]);
cP(2) = setParaFields(a{1},'Number',10,'rBin','Radial Bin Count',[2 40 true]);
cP(3) = setParaFields(a{1},'List',{1,pList},'cFunc','Expected Location Method');

% sets the tool-tip strings
cP(1).TTstr = 'Duration over which the population location is averaged';
cP(2).TTstr = 'The number of tube discretisation regions';
cP(3).TTstr = 'The function type used to calculate expected location';

% --- initialises the plotting parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 5;
pP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General','2 - Expected Location'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},'Boolean',0,'isZeitG','Use Zeitgeiber Time Format');
pP(2) = setParaFields(a{1},'Boolean',0,'plotErr','Plot Standard Error');
pP(3) = setParaFields(a{2},'Boolean',0,'plotMean','Plot Expected Location');
pP(4) = setParaFields(a{2},'Number',1,'lWid','Plot Line Width',[0.1 10 false],{2,2});
pP(5) = setParaFields(a{2},'Number',1,'dRate','Signal Downsampling Rate',[1 100 true],{2,2});
    
% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
nApp = length(snTot.iMov.ok);  
pF = setFormatFields(1);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],'');
pF.yLabel = setFormatFields([],'');
pF.Axis = setFormatFields([],[]);

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Title(i).String = snTot.iMov.pInfo.gName{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot,false,false);
[Type1,xDep1,xDep2] = deal(4,{'T','rBin'},{'T'});

% sets the independent output variables
oP = addXVarField(oP,'Time','T','Time');
oP = addXVarField(oP,'Bin Index','rBin','Group');

% sets the dependent output variables
oP = addYVarField(oP,'Prob.','P',[],5,xDep1,1);
oP = addYVarField(oP,'Location Prob.','Y',[],Type1,xDep1);
oP = addYVarField(oP,'Location (Mean)','Y_mn',[],Type1,xDep2);
oP = addYVarField(oP,'Location (SEM)','Y_sem',[],Type1,xDep2);

% ----------------------------------------------------------------------- %
% ---                       CALCULATION FUNCTION                      --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [plotD,ok] = calcFunc(snTot,pData,gPara,cP)

% initialises the calculation parameters (if not already initialised)
if (nargin == 3)
    % retrieves the parameter struct
    cP = retParaStruct(pData.cP,gPara);
end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
[Ttot,flyok] = deal(cell2mat(snTot.T),snTot.iMov.flyok);
[nApp,ok] = deal(length(flyok),true);

% determines the binned indices (for time length, tBin) and determines the
% bins which has at least two time points
indB = detTimeBinIndices(Ttot,cP.tBin);
[nFrm,T] = deal(length(indB),setupBinnedTimeArray(Ttot,indB));

% sets the bin string names
rBin = num2cell(1:cP.rBin)';
rStr = cellfun(@(x)(sprintf('Bin #%i',x)),rBin,'un',0);

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                                'T',T,'rBin',rStr,'P',[],...
                                'Y',[],'Y_mn',[],'Y_sem',[]);

% memory allocation                            
HistN = repmat({plotD(1).P},1,nApp);                            
Qw = repmat((1:cP.rBin)',1,nFrm)';               

% creates the waitbar figure
wStr = {'Overall Progress'};
h = ProgBar(wStr,'1D Location Heat Map Calculations');

% -------------------------------------- %
% --- HEATMAP HISTOGRAM CALCULATIONS --- %
% -------------------------------------- %

% determines the binned indices
h.Update(1,'Determining Time Bin Indices...',1/(2+nApp));

% loops through each of the apparatus calculating the polar coordinates
for i = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf('Calculating Histograms (%i of %i)',i,nApp);
    if h.Update(1,wStrNw,i/(1+nApp))
        % if the user cancelled, then exit the function
        ok = false; return
    end    
    
    % only calculate if values exist...
    if (~isempty(snTot.Px{i}))                      
        % calculates the position histogram
        [HistN{i},plotD(i).P] = calcPosHist(snTot,cP,indB,i);
        N = sum(snTot.iMov.flyok{i});
        
        % sets the plot values (reverses the order for the 2D values)
        x = sum(cell2mat(reshape...
                        (HistN{i},[1,1,numel(HistN{i})])),3,'omitnan');
        plotD(i).Y = x./repmat(sum(x,2),1,size(x,2));        
        
        % calculates the moving mean based on the function type
        Ex = sum(Qw.*plotD(i).Y,2,'omitnan');               
        switch cP.cFunc
            case ('Mean Location') 
                % case is the mean location
                plotD(i).Y_mn = Ex; 
                
            case ('Median Location') 
                % case is the median location
                A = num2cell(cellfun(@(x,y)(x*ones(1,y)),num2cell(Qw),...
                                num2cell(x),'un',0),2);
                plotD(i).Y_mn = cellfun(@(x)...
                                (median(cell2mat(x,'omitnan'))),A);
        end       

        % calculates the expected location standard error
        Ex2 = sum((Qw.^2).*plotD(i).Y,2,'omitnan');
        plotD(i).Y_sem = sqrt(Ex2-Ex.^2)./sqrt(N);        
    end
end

% closes the waitbar figure
if ~h.Update(1,'2D Histogram Calculations Complete!',1)
    h.closeProgBar()
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

% other parameters
[xi,dY] = deal(1:nApp,1);

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the output data
pData.pF = pF;

% retrieves the formatting struct
pF = retFormatStruct(pF,nApp);
[pF.Title,pF.yLabel.String] = deal(pF.Title(ind),'A');

% ------------------------------------- %
% --- PLOT PROPERTIES INTIALISATION --- %
% ------------------------------------- %
    
% retrieves the panel object handle
[hAx,hFig] = getCurrentAxesProp();

% sets the y-axis step size (based on the number of bins)
rBin = max(cellfun(@(x)(size(x,2)),field2cell(p,'Y')));
switch (rBin)
    case {2,3,4}
        yStep = 1;
    case {6,7,8,9,10}
        yStep = 2;
    case {11,12,13,14,15}
        yStep = 3;
    otherwise
        yStep = 4;
end

% sets the y-axis tick locations/label strings
YY = (yStep:yStep:rBin)';
yTick = cell2mat(cellfun(@(x)(YY+(x-1)*(rBin+1)),num2cell((1:nApp)'),'un',0));

% sets the y-strings based on whether or not
[yStr,dYtick] = deal(num2str(repmat(YY,nApp,1)),0);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% parameters
xLim = [1 max(cellfun(@(x)(size(x,1)),field2cell(p,'Y')))];
[xFill,yFill] = deal(xLim([1 2 2 1]),[0 0 1 1]);

% sets the entire heatmap image
CAll = NaN(nApp*rBin+(nApp-1)*dY,size(p(1).Y,1));
for i = 1:nApp    
    if (~isempty(p(i).Y))
        indNw = (1:rBin) + (i-1)*(rBin+dY);
        CAll(indNw,:) = p(i).Y';
    end
end

% creates the new image and sets the image properties
imagesc(CAll*100); axis normal; axis ij; hold on
colormap('jet'); 
caxis([0 100]); 

% updates the image y-tick positions/labels
set(gca,'ytick',yTick-dYtick,'yticklabel',yStr)
[pltInd,col] = deal(1:pP.dRate:size(CAll,2),'r');
if (pltInd(end) ~= length(p(i).Y_mn)); pltInd(end+1) = size(CAll,2); end

% loops through all the subplots 
for j = nApp:-1:1
    % sets the subplot and the an
    i = xi(j);
    yOfs = (i-1)*(rBin+1);
    set(gca,'linewidth',1.5,'box','on')

    % if not the 
    if (i ~= nApp)
        fill(xFill,yFill+(i*(rBin+1)-0.5),0.5*[1 1 1])
    end 
    
    % plots the values
    if (~isempty(p(i).Y))    
        % plots the tube row markers
        xNw = repmat(get(hAx,'xlim')',1,rBin-1);
        plot(xNw,repmat(1:rBin-1,2,1)+(0.5+yOfs),'w:')    
        
        % plots the moving mean values (if required)
        if (pP.plotErr && ~all(isnan(p(i).Y_mn)))      
            plotSignalSEM(p(i).Y_mn(pltInd)+yOfs,p(i).Y_sem(pltInd),pltInd','g',0.5)
        end          
        
        % plots the moving mean values (if required)        
        if (pP.plotMean)        
            plot(pltInd,p(i).Y_mn(pltInd)+yOfs,col,'linewidth',pP.lWid);      
        end          
    end
end

% removes the x-axis tick-labels
set(hAx,'xticklabel',[],'UserData',1)

% creates the colour bar
hCB = colorbar;

% set the figure font properties
formatPlotAxis(hAx,pF,1);
delete(findall(hAx,'tag','Title'));

% resets the axis position
resetAxesPos(hAx,1,1);
resetObjPos(hAx,'Height',-0.025,1);

% ----------------------------------- %
% --- APPARATUS NAME STRING SETUP --- %
% ----------------------------------- %

% y-label x-offset proportion
[H,hText] = deal(0.05,zeros(nApp,1));

% creates a dummy label and gets the labels position
hYLbl = get(hAx,'ylabel');
set(hYLbl,'Units','Data');
yLblP = get(hYLbl,'extent');
set(get(hAx,'ylabel'),'string',' ')

% sets the axis titles for all apparatus
for i = 1:nApp
    % creates a text object and updates the properties
    yOfs = i*(rBin+1) - 1;
    hText(i) = text(0,0,pF.Title(i).String);
    updateFontProps(hText(i),pF.yLabel.Font,i,pF.Title(i).String);
    
    % resets the location of the text object
    set(hText(i),'Rotation',90,'tag','Title','UserData',i,'Parent',hAx)
    posNw = [yLblP(1)+yLblP(3)/2 (yOfs-(rBin/2)) 0];
    set(hText(i),'position',posNw,'HorizontalAlignment','center');
end

% optimises the y-axis title placement
optYTitlePlacement(hAx,hText)

% updates the position of the colour bar
set(hCB,'Units','Normalized')
cbPos = get(hCB,'Position');
cbPos(3) = cbPos(3)*2.5;
    
% resets the axes width
resetObjPos(hAx,'width',0.99-sum(cbPos([1 3])),1)

% ---------------------------- %
% --- DAY/NIGHT AXIS SETUP --- %
% ---------------------------- %

% sets the axis height and offset
Tmlt = convertTime(1,'sec','hrs');

% retrieves the axis position (ensures units are in pixels)
if (strcmp(get(hAx,'Units'),'pixels'))
    set(hAx,'Units','Normalized')
    axPos = get(hAx,'position');
    set(hAx,'Units','Pixels')
else
    axPos = get(hAx,'position');    
end

% retrieves the current axis position
axPos2 = [axPos(1),axPos(2)+H,axPos(3),H];

% resets the axis height
resetObjPos(hAx,'height',-2*H,1)
resetObjPos(hAx,'bottom',2*H,1)
pause(0.1);

% creates the second axis for the day/night plot
set(0,'CurrentFigure',hFig);  
hAx2 = axes('parent',get(hAx,'parent'),'position',axPos2,...
            'ytick',[],'yticklabel',[],'xcolor','k','xticklabel',[],...
            'ycolor','k','linewidth',1.5,'zcolor','w','xtick',[],...
            'Units','Normalized','tag','hAxDN');
        
% plots the day/night graph onto the axis
set(hAx2,'xlim',[1 length(p(1).T)],'ylim',[0 1]+0.5)
plotDayNightAxes(hAx,snTot,p(1).T*Tmlt,1,1)

% case is using the absolute time axis
if (pP.isZeitG)            
    % case is using Zeitgeiber time
    setZeitGTimeAxis(hAx2,p(1).T,snTot);        
else
    % if plotting day/night, then set absolute time
    setAbsTimeAxis(hAx2,p(1).T,snTot);        
end

% plots the seperation markers
[xLim1,xLim2] = deal(get(hAx,'xlim'),get(hAx2,'xlim'));
[xTick2,yLim] = deal(get(hAx2,'xtick'),get(hAx,'yLim'));
xTick = (xTick2-xLim1(1))*(diff(xLim2)/diff(xLim1));
for i = 1:length(xTick)
    plot(hAx,xTick(i)*[1 1],yLim,'g-','linewidth',2)
end

% formats the axis and makes it visible
pF.Title(1).String = '';
formatPlotAxis(hAx2,pF,1);
setObjVisibility(hAx2,'on')

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)
    
% removes the flag for the population variable
pData.oP.yVar(1).Type(:) = false;

% ----------------------------------------------------------------------- %
% ---                        OTHER FUNCTIONS                          --- %
% ----------------------------------------------------------------------- %    

% --- calculates the Day/Night position histogram values --- %
function [BHistT,BHistP] = calcPosHist(snTot,cP,indB,iApp)

% parameters
[pDel,xiH] = deal(0.001,0.5:(cP.rBin+0.5));
flyok = snTot.iMov.flyok{iApp};
ifok = find(flyok);

% sets the x-coordinates of the flies   
Px = snTot.Px{iApp}(:,flyok);

% calculates the flies binned x-locations (denoted by PxN)
PxMin = repmat(min(Px,[],1)*(1-pDel),size(Px,1),1);
PxMax = repmat(max(Px,[],1)*(1+pDel),size(Px,1),1);
PxN = min(max(1,floor(cP.rBin*(Px - PxMin)./(PxMax - PxMin)) + 1),cP.rBin);
clear Px

% memory allocation
HH = NaN(length(indB),cP.rBin);
[H,BHistT] = deal(cell(length(indB),1),repmat({HH},1,length(flyok)));

% calculates the histograms for each time bin over all flies
for i = 1:length(indB)  
    if isempty(indB{i})
        H{i} = NaN(length(ifok),cP.rBin);
    else    
        PxNC = num2cell(PxN(indB{i},:),1);
        H{i} = cell2mat(cellfun(@(x)(histcounts(x,xiH)'),PxNC,'un',0))'; 
    end
end

% for each of the flies, set the histogram
for k = 1:size(PxN,2)
    j = ifok(k);
    BHistT{j} = combineNumericCells(cellfun(@(x)(x(k,:)),H,'un',0))';
end

% converts the histogram counts to percentages
BHistP = cellfun(@(x)(x./repmat(sum(x,2),1,size(x,2))),BHistT,'un',0);   
