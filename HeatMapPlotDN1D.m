% --- plots the heatmap of the daily 1D fly location over the duration of 
%     an experiment (long 1D experiment only)
function pData = HeatMapPlotDN1D(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = '1D Fly Location Heatmap (Long)';
pData.Type = {'Pop','Multi'};
pData.fType = [2 1 3 2];
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
rI.Dur = 'Long';
rI.Shape = '1D';
rI.Stim = 'None';
rI.Spec = 'None';
        
% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
Tmin = 3;                       % short experiment max duration (in hours)
nPara = 4;
cP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the parameter list for 
pList = {'Mean Location','Median Location'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'Number',60,'tBin','Time Bin Size (sec)',[15 600 false]);
cP(2) = setParaFields(a{1},'Number',10,'rBin','Region Discretisation Count',[2 40 true]);
cP(3) = setParaFields(a{1},'List',{1,pList},'cFunc','Expected Location Method');

% sets the tool-tip strings
cP(1).TTstr = 'Duration over which the population location is averaged';
cP(2).TTstr = 'The number of tube discretisation regions';
cP(3).TTstr = 'The function type used to calculate expected location';

% if short experiments only, then create alignment checkbox option
if (all(convertTime(cellfun(@(x)(x{end}(end)),field2cell(snTot,'T')),'sec','hrs') <= Tmin))
    % creates option
    cP(4) = setParaFields([],'Boolean',0,'isAlign','Align All Experiments Start Points');
    cP(4).TTstr = 'Aligns all experiments to their respective start times (Short experiments only)';    
else
    % long experiment, so remove option
    cP = cP(1:3);
end

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
pP(4) = setParaFields(a{2},'Number',1,'lWid','Plot Line Width',[0.1 10 false],{3,2});
pP(5) = setParaFields(a{2},'Number',1,'dRate','Signal Downsampling Rate',[1 100 true],{3,2});
    
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
oP = setupOutputParaStruct(snTot);
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
[nApp,ok] = deal(length(snTot(1).iMov.ok),true);
nExp = length(snTot);

% sets the alignment flag
if (isfield(cP,'isAlign'))
    isAlign = cP.isAlign;
else
    isAlign = false;
end

% determines if signals are being aligned or not...
Tf = convertTime(cellfun(@(x)(x{end}(end)),field2cell(snTot,'T')),'sec','hours');
if (~isAlign)
    % retrieves the final times for each of the experiments (in days)     
    Tmax = convertTime(1,'day','sec');  
else
    Tmax = convertTime(max(Tf),'hours','sec');    
end
T = ((cP.tBin/2):(cP.tBin):Tmax)';

% sets the bin string names
rBin = num2cell(1:cP.rBin)';
bStr = cellfun(@(x)(sprintf('Bin #%i',x)),rBin,'un',0);

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                                'T',T,'rBin',bStr,'P',[],...
                                'Y',[],'Y_mn',[],'Y_sem',[]);

% creates the waitbar figure
wStr = {'Overall Progress','Location Heatmap Calculations'};
wOfs = (nExp > 1); wStr = wStr((2-wOfs):end);
h = ProgBar(wStr,'2D Position Heatmap Calculations');

% fixed parameters
tDayStart = convertTime(cP.Tgrp0,'hours','secs');     % sets the day start offset
nDay = max(detExptDayDuration(snTot));

% memory allocation
nFrm = length(T);
HistN = repmat({plotD(1).P},1,nApp);
N = zeros(nFrm,nApp);
Qw = repmat((1:cP.rBin)',1,nFrm)';     

% -------------------------------------- %
% --- HEATMAP HISTOGRAM CALCULATIONS --- %
% -------------------------------------- %

% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure (if more than one solution file)
    if wOfs > 0
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
        if h.Update(1,wStrNw,i/(1+nExp))
            [plotD,ok] = deal([],false);
            return
        else
            h.Update(2,wStr{2},0);
        end
    end    
    
    % sets the experiment start times
    if isAlign
        T0nw = [0,cP.Tgrp0,0,0];
    else
        T0nw = [0,snTot(i).iExpt(1).Timing.T0(4:end)];
    end                
    
    % determines the index offset    
    Ttot = cell2mat(snTot(i).T);             

    % determines the binned indices
    h.Update(1+wOfs,'Determining Time Bin Indices...',1/(2+nApp));        
    indB = detTimeBinIndices(Ttot,cP.tBin); 

%     %
%     if (length(indB)) > (length(T))
%         indB = indB(1:length(T));
%     end
    
    % sets averaging indices for the data
    nFrmMx = min(length(indB),ceil(Ttot(end)/cP.tBin));
    indOfs = floor(mod((vec2sec(T0nw)-tDayStart)/cP.tBin,nFrm));
    ii = arrayfun(@(x)((x-indOfs):nFrm:nFrmMx),(1:nFrm)','un',0);
    ii = cellfun(@(x)(x(x>0)),ii,'un',0);   
        
    % loops through each apparatus plotting the position heat-maps
    for j = 1:nApp
        % if the user cancelled, then exit the function
        wStrNw = sprintf('Calculating Histograms (%i of %i)',j,nApp);
        if h.Update(1+wOfs,wStrNw,j/(1+nApp))
            return
        end
                
        % only calculate if values exist...
        if ~isempty(snTot(i).Px)
            [HistNnw,HistPnw] = calcDNPosHist(snTot(i),nDay,cP,indB,ii,j); 
            if ~isempty(HistNnw)
                jj = 1:size(HistNnw,2);
                HistN{j}(:,jj,i) = HistNnw;
                plotD(j).P(:,jj,i) = HistPnw;
                N(:,j) = N(:,j) + sum(snTot(i).iMov.flyok{j});
            end
        end
    end
end

% calculates the final values for all
for j = 1:nApp
    % sets the plot values (reverses the order for the 2D values)
    x = sum(cell2mat(reshape(HistN{j},[1,1,numel(HistN{j})])),3,'omitnan');
    plotD(j).Y = x./repmat(sum(x,2),1,size(x,2));

    % calculates the moving mean based on the function type
    Ex = sum(Qw.*plotD(j).Y,2,'omitnan');
    switch cP.cFunc
        case ('Mean Location') % case is the mean
            plotD(j).Y_mn = Ex;                         
        case ('Median Location') % case is the median
            A = num2cell(cellfun(@(x,y)(x*ones(1,y)),num2cell(Qw),...
                                num2cell(x),'un',0),2);
            plotD(j).Y_mn = cellfun(@(x)...
                                (median(cell2mat(x),'omitnan')),A);
    end      

    % calculates the expected location standard error
    Ex2 = sum((Qw.^2).*plotD(j).Y,2,'omitnan');
    plotD(j).Y_sem = sqrt(Ex2-Ex.^2)./sqrt(N(:,j));
end
        
% closes the waitbar figure
if ~h.Update(1,'Histogram Calculations Complete!',1)
    h.closeProgBar()
end

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

% sets the plotting indices and subplot indices
ind = find(sP.Sub.isPlot);
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% memory allocation
[xi,dY] = deal(1:nApp,1);

% ------------------------------------- %
% --- PLOT PROPERTIES INTIALISATION --- %
% ------------------------------------- %
    
% sets the axis handle
[hAx,hFig] = getCurrentAxesProp;

% sets the y-axis step size (based on the number of bins)
rBin = size(p(1).Y,2);
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

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the output data
pData.pF = pF;

% retrieves the formatting struct
pF = retFormatStruct(pF,nApp);
[pF.Title,pF.yLabel.String] = deal(pF.Title(ind),'A');

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% sets the entire heatmap image
CAll = NaN(nApp*rBin+(nApp-1)*dY,size(p(1).Y,1));
for i = 1:nApp    
    indNw = (1:rBin) + (i-1)*(rBin+dY);
    CAll(indNw,:) = p(i).Y';
end

% creates the new image and sets the image properties
imagesc(100*CAll); axis normal; axis ij; hold on
colormap('jet'); 
caxis([0 100]); 

% parameters
xLim = getCurrentAxesProp('xlim');
[xFill,yFill] = deal(xLim([1 2 2 1]),[0 0 1 1]);

% updates the image y-tick positions/labels
set(gca,'ytick',yTick-dYtick,'yticklabel',yStr)
[pltInd,col] = deal(1:pP.dRate:length(p(i).Y_mn),'r');
if (pltInd(end) ~= length(p(i).Y_mn)); pltInd(end+1) = length(p(i).Y_mn); end

% loops through all the subplots 
for j = 1:nApp
    % sets the subplot and the an
    i = xi(j);
    yOfs = (i-1)*(rBin+1);
    set(gca,'linewidth',1.5,'box','on')

    % if not the final heat-map, then create a seperator bar
    if (i ~= nApp)
        fill(xFill,yFill+(i*(rBin+1)-0.5),0.5*[1 1 1])
    end
        
    % plots the tube row markers    
    plot(repmat(get(hAx,'xlim')',1,rBin-1),repmat(1:rBin-1,2,1)+(0.5+yOfs),'w:')    

    % plots the moving mean values (if required)
    if (pP.plotErr && ~all(isnan(p(i).Y_mn)))      
        plotSignalSEM(p(i).Y_mn(pltInd)+yOfs,p(i).Y_sem(pltInd),pltInd','g',0.5)
    end              
    
    % plots the moving mean values (if required)
    if (pP.plotMean)        
        plot(pltInd,p(i).Y_mn(pltInd)+yOfs,col,'linewidth',pP.lWid);      
    end          
end

% removes the x-axis tick-labels
set(hAx,'xticklabel',[],'UserData',1,'ticklength',[0 0])

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

% sets the day/night axis markers
[ix,iy,yLim,fAlpha,dT] = deal([1 1 2 2],[1 2 2 1],0.5+[0 1],0.9,0.01);
[xFillD,xFillN,yFill] = deal([0 12],[12 24],yLim);

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
hAx2 = axes('parent',get(hAx,'parent'),'Units','Normalized',...
            'ytick',[],'yticklabel',[],'xcolor','k','xticklabel',[],...
            'ycolor','k','linewidth',1.5,'zcolor','w','xtick',[],...
            'visible','off','position',axPos2);
set(hFig,'CurrentAxes',hAx2);
pause(0.05);

% adds in the day/night fill objects (first trace only)
hold(hAx2,'on'); axis(hAx2,'on')
hFD = fill(xFillD(ix),yFill(iy),0.9*[1 1 0],'tag','hDN');
hFN = fill(xFillN(ix),yFill(iy),0.5*[1 1 1],'tag','hDN');  

set(hFD,'parent',hAx2); 
set(hFN,'parent',hAx2)

% sets the axis limits
set(hAx2,'xlim',[0 24],'ylim',yLim);
set(hAx,'xlim',convertTime([0 24],'hours','min'))

% case is using the absolute time axis
Tax = convertTime([0 24],'hrs','sec');
if (pP.isZeitG)            
    % case is using Zeitgeiber time
    setZeitGTimeAxis(hAx2,Tax);        
else
    % if plotting day/night, then set absolute time
    setAbsTimeAxis(hAx2,Tax);        
end    

% formats the axis and makes it visible
pF.Title(1).String = '';
formatPlotAxis(hAx2,pF,1);
expandAxesLim(hAx2);
    
% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% % removes the flag for the population variable
% pData.oP.yVar(1).Type(:) = false;

% ----------------------------------------------------------------------- %
% ---                        OTHER FUNCTIONS                          --- %
% ----------------------------------------------------------------------- %    

% --- calculates the Day/Night position histogram values --- %
function [BHistT,BHistP] = calcDNPosHist(snTot,nDay,cP,indB,ii,iApp)

% parameters
[nBin,nGrp] = deal(convertTime(1,'day','sec')/cP.tBin, length(indB));
[pDel,xiH,flyok] = deal(0.001,0.5:(cP.rBin+0.5),snTot.iMov.flyok{iApp});
ifok = find(flyok);

%
if (isempty(ifok))
    [BHistT,BHistP] = deal([]);
    return
end

% % determines the 
% jj = find(~cellfun('isempty',ii));
% ii = ii(jj);
% 
% % sets the day indices 
% [iR,iC] = deal(mod(jj-1,nBin)+1,floor((jj-1)/nBin)+1);
% indT = NaN(nBin,max(iC));
% indT(sub2ind(size(indT),iR,iC)) = 1:length(iR);

% determines the 
jj = find(~cellfun('isempty',ii));

% inserts NaN's to make sure the day's align correctly
i0 = jj(cellfun(@(x)(x(1)),ii(jj))==1);

%
kk = 1:(i0-1);
isNE = kk(~cellfun('isempty',ii(kk)));
ii(isNE) = cellfun(@(x)([NaN,x]),ii(isNE),'un',0);

% converts the indices into 
indT = combineNumericCells(ii)';
if (size(indT,1) < nBin)
    indT = [indT;NaN(nBin-size(indT,1),size(indT,2))];
end
 
% ensures there are no zeros after the NaN's in the last column
isNN_L = find(isnan(indT(:,end)),1,'last');
if ~isempty(isNN_L)
    if isNN_L ~= size(indT, 1)
        indT(isNN_L:end, end) = NaN;
    end
end
    
% memory allocation
nDay = max(nDay,size(indT,2));
[HH,indT] = deal(NaN(nBin,cP.rBin),num2cell(indT,1)');
[H,BHistT] = deal(cell(nGrp,1),repmat({HH},nDay,length(flyok)));

% sets the x-coordinates of the flies   
Px = snTot.Px{iApp}(:,flyok);

% calculates the flies binned x-locations (denoted by PxN)
PxMin = repmat(min(Px,[],1)*(1-pDel),size(Px,1),1);
PxMax = repmat(max(Px,[],1)*(1+pDel),size(Px,1),1);
PxN = min(max(1,floor(cP.rBin*(Px - PxMin)./(PxMax - PxMin)) + 1),cP.rBin);
clear Px;

% calculates the histograms for each time bin over all flies
for i = 1:nGrp
    if (i > length(indB)) || isempty(indB{i})
        H{i} = NaN(length(ifok),cP.rBin);
    else
        PxNC = num2cell(PxN(indB{i},:),1);
        H{i} = cell2mat(cellfun(@(x)(histcounts(x,xiH)'),PxNC,'un',0))';         
    end
end

% for each of the flies, set the histogram
for k = 1:size(PxN,2) 
    [j,ii] = deal(ifok(k),1:length(indT));
    Htmp = combineNumericCells(cellfun(@(x)(x(k,:)),H,'un',0))';
    BHistT(ii,j) = cellfun(@(x)(setHistArray(Htmp,x)),indT,'un',0);    
end 

% converts the histogram counts to percentages
BHistP = cellfun(@(x)(x./repmat(sum(x,2),1,size(x,2))),BHistT,'un',0);

% --- sets the histogram values for each time group
function BHist = setHistArray(H,ind)

% memory allocation
[BHist,isNN] = deal(NaN(length(ind),size(H,2)),~isnan(ind));
ii = find(isNN,1,'first'):find(isNN,1,'last');

% sets the histogram values
BHist(ii,:) = H(ind(ii),:);
% BHist(ii,:) = H(ii,:);
