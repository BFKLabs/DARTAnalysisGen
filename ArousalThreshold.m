% --- calculates the level of fly arousability with respect to delivered 
%     stimuli throughout the day (long experiment only)
function pData = ArousalThreshold(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);

% sets the function name/type
pData.Name = 'Arousal Threshold';
pData.Type = {'Pop','Multi'};
pData.fType = [3 2 3 1 0];

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
    [pData.hasSP,pData.hasRS] = deal(true,false);
    pData.sP = initSpecialPara(snTotL,pData);
end

% ----------------------------------------------------------------------- %
% ---                 PARAMETER STRUCT SETUP FUNCTIONS                --- %
% ----------------------------------------------------------------------- %

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 1;
cP = setParaFields(nPara);
isLong = ~detIfShortExpt(field2cell(snTot,'T'));

% sets the tab list names
a = {'1 - General'};

% sets the parameter list for
pList = cellfun(@num2str,num2cell([1 2 4 6 8 12 24]),'un',0);

% sets the parameter fields
cP(1) = setParaFields(a{1},'List',{1+6*isLong,pList},...
                    'nGrp','Number of Daily Time Groups',[],{0,~isLong});

% sets the tool-tip strings
cP(1).TTstr = 'The number of groups that the day is split up into';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 6;
pP = setParaFields(nPara);

% sets the parameter list for 
pList = {'Box & Whisker','CDF Heatmap'};

% sets the tab list names
a = {'1 - General','2 - Quartile Levels'};

% sets the plot parameter fields into the data struct
pP(1) = setParaFields(a{1},'List',{1,pList},'pType','Arousal Threshold Plot Type');
pP(2) = setParaFields(a{1},'Number',2,'lWid','Plot Line Width',[0.5 10 false],{1,2});
pP(3) = setParaFields(a{1},'Boolean',0,'incMoving','Include Pre-Stimuli Moving Flies');
pP(4) = setParaFields(a{1},'Boolean',1,'pltDN','Plot Day/Night Background Image',[],{1,1});
pP(5) = setParaFields(a{2},'Boolean',0,'pltMedian','Plot Median Level',[],{1,2});
pP(6) = setParaFields(a{2},'Boolean',0,'pltQuart','Plot Inter-Quartile Levels',[],{1,2});

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% determines the apparatus count
nApp = length(snTot.appPara.ok);  
pF = setFormatFields(nApp);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],'');
pF.yLabel = setFormatFields([],'',1);
pF.Axis = setFormatFields([],[]);

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Title(i).String = snTot.appPara.Name{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);
[Stats,Type1,Type2] = deal({'CompMulti','Tgrp'},5,4);
[xDep2,xDep3,xDep4] = deal({'Tgrp'},{'Tgrp','yStr'},{'Tgrp','yStrIm'});

% sets the independent output variables
oP = addXVarField(oP,'Time Group','Tgrp','Other');
oP = addXVarField(oP,'Stimuli Count','iStim','Other');
oP = addXVarField(oP,'React Lvl (All)','yStr','Other');
oP = addXVarField(oP,'React Lvl (Immob)','yStrIm','Other');

% sets the dependent output variables
oP = addYVarField(oP,'Arousal Lvl','indR',Stats,Type1,xDep2,1);
oP = addYVarField(oP,'Arousal Lvl (Med)','indR_md',[],Type2,xDep2);
oP = addYVarField(oP,'Arousal Lvl (LQ)','indR_lq',[],Type2,xDep2);
oP = addYVarField(oP,'Arousal Lvl (UQ)','indR_uq',[],Type2,xDep2);
oP = addYVarField(oP,'Immob Arousal Lvl','indRIm',Stats,[],xDep2,1);
oP = addYVarField(oP,'Immob Arousal Lvl (Med)','indRIm_md',[],Type2,xDep2);
oP = addYVarField(oP,'Immob Arousal Lvl (LQ)','indRIm_lq',[],Type2,xDep2);
oP = addYVarField(oP,'Immob Arousal Lvl (UQ)','indRIm_uq',[],Type2,xDep2);
oP = addYVarField(oP,'Lvl CDF (All)','pCDF',[],Type2,xDep3,1);
oP = addYVarField(oP,'Lvl CDF (Immob)','pCDFim',[],Type2,xDep4,1);

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

% determines the currently loaded experimental file is valid for
% calculating arousal threshold
if (length(snTot(1).iExpt.Stim(1).sigPara{1}) == 1)
    eStr = {'Error! Single-phase only stimuli events detected within data set.';
        'To calculate arousal threshold, choose a multi-phase stimuli experiment.'};
    waitfor(errordlg(eStr,'Not A Multi-Phase Experiment')); 
    [plotD,ok] = deal([],false); return
end

% sets the movement calculation type
cP.movType = 'Absolute Speed';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
[nApp,nExp,ok] = deal(length(snTot(1).appPara.flyok),length(snTot),true);
nGrp = str2double(cP.nGrp);

% fixed parameters
a = NaN(nGrp,1);
T = 1:nGrp;

% sets the daily time group strings
Tgrp = setTimeGroupStrings(nGrp,cP.Tgrp0);  
iStim = cellfun(@(x)(sprintf('Stim #%i',x)),num2cell(1:24/length(Tgrp))','un',0);

% sets the y-axis level strings
Nlvl = length(snTot(1).iExpt.Stim(1).sigPara{1})+2;
yStr = [{'MF'};cellfun(@(x)(sprintf('%.1f',x)),num2cell(...
                linspace(0,1,Nlvl-2))','un',0);{'NR'}];                                  

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                                 'Tgrp',Tgrp,'pCDF',a,'pCDFim',a,'T',T,...
                                 'yStr',yStr,'yStrIm',yStr(2:end),'iStim',iStim,...
                                 'indR',[],'indRIm',[],'indRT',[],'indRImT',[],...
                                 'indR_md',[],'indR_lq',[],'indR_uq',[],...
                                 'indRIm_md',[],'indRIm_lq',[],'indRIm_uq',[]);

% ---------------------------------------------------- %
% --- STIMULI TRACE & IMMOBILITY TIME CALCULATIONS --- %
% ---------------------------------------------------- %

% sets the waitbar offset (is > 0 for more than one
wStr = {'Overall Progress','Calculating Arousal Threshold Levels'};
wOfs = (length(snTot) > 1);

% creates the waitbar figure
wStr = wStr((2-wOfs):end);
h = ProgBar(wStr,'Arousal Threshold Calculations');

% loops through each of the
for i = 1:nExp
    % updates the waitbar figure (if more than one solution file)
    if wOfs > 0
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
        h.Update(1,wStrNw,i/(1+nExp));
    end
    
    % sets the video/stimuli time stamps into a single vector
    if ~isempty(cell2mat(snTot(i).Ts))        
        % calculates the new sleep intensity data
        %     [YcountNw,YcountRNw,~,TreactNw] = ...
        indReactNw = calcArousalThresholdData(snTot(i),cP,h,wOfs);
        if (isempty(indReactNw))
            % if the user cancelled, then exit the function
            [plotD,ok] = deal([],false);
            return
        else
            % otherwise, append the data to the arrays
            for j = 1:nApp
                if (snTot(i).appPara.ok(j))
                    % sets the rows/column indices
                    xiR = 1:size(indReactNw{j},1);
                    xiC = 1:size(indReactNw{j},2);
                    
                    % sets the raw data for all flies
                    plotD(j).indR(xiR,xiC,i) = indReactNw{j};
                    
                    % sets all the pre-stim moving fly to NaN values
                    for k = 1:numel(indReactNw{j})
                        isZ = indReactNw{j}{k} == 0;
                        indReactNw{j}{k}(isZ) = NaN;
                    end
                    
                    % sets the raw data for the immobile flies
                    plotD(j).indRIm(xiR,xiC,i) = indReactNw{j};                    
                end
            end
        end
    end
end

% calculates the total/immobile arousal threshold CDF distributions over
% all the apparatus
for i = 1:nApp    
    % calculates the metric statistics
    plotD(i) = calcMetricStats(plotD(i),'indR',4);
    plotD(i) = calcMetricStats(plotD(i),'indRIm',4);
end
 
% removes any rows/columns which are empty over all groups
plotD = removeEmptyElements(plotD, {'indR', 'indRIm'});

%
for i = 1:nApp
    % sets the combined raw data arrays
    plotD(i).indRT = cell2mat(plotD(i).indR(:));
    plotD(i).indRImT = cell2mat(plotD(i).indRIm(:));

    % calculates the empirical cdfs for each time group
    plotD(i).pCDF = cell2mat(cellfun(@(x)(calcThresholdCDF(...
                x,Nlvl)),num2cell(plotD(i).indRT,1),'un',0));
    plotD(i).pCDFim = cell2mat(cellfun(@(x)(calcThresholdCDF(...
                x,Nlvl-1,1)),num2cell(plotD(i).indRImT,1),'un',0));                                               
end

% closes the waitbar figure
if ~h.Update(1,'Arousal Threshold Calculations Complete!',1)
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
[ind,m,n] = deal(find(sP.Sub.isPlot),sP.Sub.nRow,sP.Sub.nCol);
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% sets the time group number
nGrp = str2double(cP.nGrp);

% other parameters
fAlpha = 0.4;
[ix,iy] = deal([1 1 2 2],[1 2 2 1]);

% calculates the upper y-axis limits
nLvl = size(p(1).pCDF,1) - 2;

% retrieves the panel object handle
hP = get(gca,'Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% retrieves the formatting struct
if (isempty(m)); szMx = 1; else; szMx = max([m n]); end
pF = retFormatStruct(pF,szMx);

% sets the y-axis strings
yStr = [cellfun(@(x)(sprintf('%.1f',x)),num2cell(...
                        linspace(0,1,nLvl))','un',0);{'NR'}];       
if (pP.incMoving); yStr = [{'MF'};yStr]; end            
            
% sets the y-axis label string
pF.yLabel.String = 'Arousal Threshold (g)'; 
[pF.xLabel.ind,pF.yLabel.ind] = deal(NaN);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% memory allocation
hAx = cell(nApp,1);

% day/night background fill object coordinates
xLim = [1 nGrp]+0.5*[-1 1];
[xFillD,xFillN,yFill] = deal([xLim(1) mean(xLim)],[mean(xLim) xLim(2)],[-1 nLvl+2]);

% loops through all the indices plotting the solutions
for i = 1:nApp
    % initialises the subplot axes
    hAx{i} = createSubPlotAxes(hP,[m,n],i);
    
    % creates the plot based on the plot type
    switch (pP.pType)
        case ('Box & Whisker') % case is box and whisker plot
            % plots the background image
            if (pP.pltDN && (nGrp > 1))       
                fill(xFillD(ix),yFill(iy),'y','FaceAlpha',fAlpha','tag','hDN')
                fill(xFillN(ix),yFill(iy),'k','FaceAlpha',fAlpha','tag','hDN')        
            end               
            
            % sets the box/whisker values
            if (~isempty(p(i).indRT))
                if (pP.incMoving)
                    Z = p(i).indRT; 
                else
                    Z = p(i).indRImT;
                end

                % creates the boxplot and removes the outliers
                h = boxplot(hAx{i},Z,'labels',repmat({''},1,size(Z,2)),'boxstyle','filled'); 
                set(findall(h,'tag','Median'),'Linewidth',2)
                delete(findall(h,'tag','Outliers'))            
            end
            
            % sets the axis properties
            iOfs = (1-pP.incMoving);
            set(hAx{i},'ytick',(0:length(yStr)-1)'+iOfs,'yticklabel',yStr,...
                       'ylim',[0 length(yStr)]-(0.5-iOfs));
            set(hAx{i},'ygrid','on')
        case ('CDF Heatmap') % case is the CDF heatmaps                
            % sets the heatmap values
            if (pP.incMoving) 
                % case is the including all flies 
                [Z,Ymd] = deal(p(i).pCDF,p(i).indR_md);
                [Ylu,Yuq] = deal(p(i).indR_lq,p(i).indR_uq);
            else
                % case is the including only the pre-stimuli immobile flies 
                [Z,Ymd] = deal(p(i).pCDFim,p(i).indRIm_md);
                [Ylu,Yuq] = deal(p(i).indRIm_lq,p(i).indRIm_uq);
            end
            
            % creates the image and sets the colourmap/colourbar
            imagesc(Z*100); 
            set(hAx{i},'clim',[0 100])
            colormap('jet'); 
                        
            % resets the axis limits            
            set(hAx{i},'xlim',[0 nGrp]+0.5,'ylim',[0 size(Z,1)]+0.5,...
                        'ytick',(1:length(yStr))','yticklabel',yStr)
                    
            % plots the mean (if checked)
            xi = 0:(nGrp+1);
            if (pP.pltMedian)
                % wraps around the median values
                Ymd = max(0,[mean(Ymd([1 end])),Ymd,mean(Ymd([1 end]))]);
                
                % plots the median values
                plot(hAx{i},xi,Ymd,'k','linewidth',pP.lWid)
            end
            
            % plots the inter-quartile range (if checked)
            if (pP.pltQuart)
                % wraps the lower/upper quartile values
                Ylu = max(0,[mean(Ylu([1 end])),Ylu,mean(Ylu([1 end]))]);
                Yuq = max(0,[mean(Yuq([1 end])),Yuq,mean(Yuq([1 end]))]);
                
                % plots the lower/upper quartile values
                plot(hAx{i},xi,Ylu,'k:','linewidth',pP.lWid)                
                plot(hAx{i},xi,Yuq,'k:','linewidth',pP.lWid)                
            end
    end

    % sets the time markers (if there are enough groups)
    if (nGrp > 1)                
        % sets the time markers
        NN = 2*(1 + (nGrp >= 4));
        xL = get(hAx{i},'xlim');
        dT = diff(xL); [T,T0] = deal(xL(1) + dT*(0:NN)/NN); 

        % sets the x-axis tick labels
        T = (24/nGrp)*(T-0.5)+0.5;
        T([1 end]) = T([1 end]) + 0.01*[1 -1];                                
        setAbsTimeAxis(hAx{i},convertTime(T,'hrs','sec'));                 

        % resets the tick-marks
        set(hAx{i},'xtick',T0)
    end    
    
    % sets the axis formats
    formatPlotAxis(hAx{i},pF,ind(i),1);         
end

% formats and resets the axis positions
formatMultiXYLabels(hAx,pF,[m,n]);
resetAxesPos(hAx,m,n,[0.01 0]);

% creates a colour bar for the CDF heatmaps
if (strcmp(pP.pType,'CDF Heatmap'))
    % creates the colour bar
    hCB = colorbar(hAx{end});    
    set(hCB,'Units','Normalized');

    % updates the position of the colour bar    
    if (m > 1)
        % retrieves the width of the colour bar axis
        if (isHG1)
            cbPos = get(hCB,'OuterPosition');    
            wCB = cbPos(3);
        else
            cbPos = get(hCB,'Position');
            hExt = cell2mat(cellfun(@(x)(getTextSize(hAx{1},x)),hCB.TickLabels,'un',0));
            wCB = cbPos(3) + 0.75*max(hExt(:,3));
        end    
        
        % resets the location of the final axes and colour bar position
        resetObjPos(hAx{end},'width',retObjDimPos(hAx(1),3,1));
        resetLegendAxisPos(hAx,hCB,[m n],0.005,wCB)
    end    
end

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% resets the data arrays for each group
for i = 1:length(plotD)
    plotD(i).indR = cellfun(@(x)(x'),plotD(i).indR,'un',0);
    plotD(i).indRIm = cellfun(@(x)(x'),plotD(i).indRIm,'un',0);    
    [plotD(i).pCDF,plotD(i).pCDFim] = deal(plotD(i).pCDF',plotD(i).pCDFim');
end

% resets the day separation flag (if there are no expts > 1 day in length)
pData.oP.sepDay = any(detExptDayDuration(snTot) > 1);

% ----------------------------------------------------------------------- %
% ---                         OTHER FUNCTIONS                         --- %
% ----------------------------------------------------------------------- %

% --- calculates the significant movement arousal threshold level
function Pstr = getSignifLevelStr(P)

% ensures the proportions/sample size arrays are vectors
pLvl = [log10(0.05) -2 -3 -5 -10 -20 -50 -inf];
ii = find(log10(P) <= pLvl,1,'last'); 
if (isempty(ii)); ii = 0; end

% sets the p-level signficance strings
switch (ii)
    case(0)
        Pstr = {'N/S'};
    case(1)
        Pstr = {'0.05'}; 
    case(2)
        Pstr = {'0.01'}; 
    case(3)
        Pstr = {'0.001'}; 
    otherwise
        Pstr = {sprintf('10^(%i)',pLvl(ii))}; 
end

% --- calculates the cdf proportions for N levels --- %
function pCDF = calcThresholdCDF(y,N,yMin)

% sets the offset to zero if not provided
if (nargin == 2); yMin = 0; end

% memory allocation
pCDF = ones(N,1);
if (all(y <= yMin)) || (all(isnan(y))); return; end

% calculates the ecdf and set the probabilities
[f,x] = ecdf(y(y>=yMin)); 

% sets the final array
pCDF(1:(x(1)-1)) = 0;
for i = 1:(length(x)-1)
    pCDF((x(i):(x(i+1)-1))+(1-yMin)) = f(i);
end
