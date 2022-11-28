% --- calculates the fly activity initiation (number of starts per second 
%     stopped) over the duration of an experiment (short experiment only)
function pData = IBIMetrics(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Inter-Bout Interval Statistics';
pData.Type = {'Pop','Multi'};
pData.fType = [2 1 1 1];
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
rI.Shape = 'None';
rI.Stim = 'None';
rI.Spec = 'None';
        
% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 2;                       
cP = setParaFields(nPara);

% sets the parameter fields
cP(1) = setParaFields([],'Number',2,'vAct','Activity Threshold (mm/s)',[0.1 10 false]);
cP(2) = setParaFields([],'Number',0.2,'tMove','Movement Start Duration (s)',[0.1 10 false]);

% sets the tool-tip strings
cP(1).TTstr = 'The inter-frame velocity threshold used to indicate movement';
cP(2).TTstr = 'The duration the fly has to move to be considered a start';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 7;
pP = setParaFields(nPara);

% sets the parameter list strings
pList = {'Mean Survival Curve','Shape Factor','Scale Factor',...
         'Burstiness','Memory','Burstiness + Memory'};
pList2 = {'Boxplot','Bar Graph'};

% parameter tab strings
a = {'1 - General','2 - Survival Curve','3 - Other Metrics'};

% sets the parameter fields
pP(1) = setParaFields(a{1},'List',{1,pList},'pMet','Plot Metric');
pP(2) = setParaFields(a{3},'List',{1,pList2},'pType','Graph Plot Type',[],{1,2:4});
pP(3) = setParaFields(a{2},'Number',1,'pL','Plot Line Thickness',[0.5 5 false],{1,1});
pP(4) = setParaFields(a{3},'Number',0.75,'pW','Bar Plot Relative Width',[0 1 false],{2,2});
pP(5) = setParaFields(a{1},'Boolean',1,'plotGrid','Show Axis Gridlines');
pP(6) = setParaFields(a{2},'Boolean',1,'plotFit','Plot Survival Curve Fitted Response',[],{1,1});
pP(7) = setParaFields(a{3},'Boolean',1,'plotErr','Show Error Bars/Outliers',[],{1,2:4});

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
pF = setFormatFields(1);

% initialises the font structs
pF.Title = setFormatFields([],'',1);
pF.xLabel = setFormatFields([],'');
pF.yLabel = setFormatFields([],'');
pF.Axis = setFormatFields([],[]);

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);
[Type1,Type2,Type3,xDep,Stats] = deal(3,2,4,{'xB'},{'Comp'});
oP.sepDay = false;

% sets the independent output variables
oP = addXVarField(oP,'Duration (sec)','xB','Other');

% sets the dependent output variables
oP = addYVarField(oP,'Survival Curve (Raw)','YmnR',[],Type3,xDep);
oP = addYVarField(oP,'Survival Curve (Fitted)','YmnF',[],Type3,xDep);
oP = addYVarField(oP,'Shape Factor','k',Stats,Type1,[],1);
oP = addYVarField(oP,'Shape Factor (Mean)','k_mn',[],Type2);
oP = addYVarField(oP,'Shape Factor (SEM)','k_sem',[],Type2);
oP = addYVarField(oP,'Scale Factor','l',Stats,Type1,[],1);
oP = addYVarField(oP,'Scale Factor (Mean)','l_mn',[],Type2);
oP = addYVarField(oP,'Scale Factor (SEM)','l_sem',[],Type2);
oP = addYVarField(oP,'Burstiness','B',Stats,Type1,[],1);
oP = addYVarField(oP,'Burstiness (Mean)','B_mn',[],Type2);
oP = addYVarField(oP,'Burstiness (SEM)','B_sem',[],Type2);
oP = addYVarField(oP,'Memory','M',Stats,Type1,[],1);
oP = addYVarField(oP,'Memory (Mean)','M_mn',[],Type2);
oP = addYVarField(oP,'Memory (SEM)','M_sem',[],Type2);
                         
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

% sets the movement calculation type
cP.movType = 'Absolute Speed';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
[nApp,nExp,ok] = deal(length(snTot(1).iMov.flyok),length(snTot),true);

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                              'xB',[],'YmnR',[],'YmnF',[],'kmn',[],'lmn',[],...
                              'k',[],'k_mn',[],'k_sem',[],...
                              'l',[],'l_mn',[],'l_sem',[],...
                              'B',[],'B_mn',[],'B_sem',[],...
                              'M',[],'M_mn',[],'M_sem',[]);

% other memory allocations
indT = cell(nExp,1);
[I,TT] = deal(cell(nApp,nExp),cell(1,nExp));

% creates the waitbar figure
wStr = {'Movement Detection','Inter-Bout Interval Calculations'};
h = ProgBar(wStr,'Activity Calculations');

% determines the maximum fly count
nFly = zeros(nApp,1);
for i = 1:length(snTot)
    % sets the number of flies in the experiment
    nFly = max(nFly,cellfun('length',snTot(i).iMov.flyok));
end

% ------------------------------------------------------- %
% --- INTER-FRAME DISTANCE CALCULATION & THRESHOLDING --- %
% ------------------------------------------------------- %
    
% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure 
    wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
    if h.Update(1,wStrNw,i/nExp)
        [plotD,ok] = deal([],false);
        return
    end

    % calculates the video frame rate and experiment apparatus indices
    FPS = 1/mean(diff(T{i}),'omitnan');    
    iApp = find(~cellfun('isempty',snTot(i).iMov.flyok));
    
    % determines the points where there is a large time gap
    Tp = snTot(i).iExpt.Timing.Tp*3;
    indT{i} = find(diff(T{i}) > Tp);
    
    % sets the relevant time points and apparatus indices for this expt
    ii = 1:length(T{i});
    isMove = calcFlyMove(snTot(i),T{i},ii,iApp,cP.vAct);      
    
    % calculates the time mid point of each frame, and removes all the
    % frames where the inter-frame time difference is large (ie, the
    % change-over in video recordings)
    [Tmid,isOK] = deal((T{i}(1:end-1) + T{i}(2:end))/2,diff(T{i}) < 2/FPS);
    
    % sets the 
    TT{i} = Tmid(isOK);
    for j = 1:length(isMove)
        I{iApp(j),i} = isMove{j}(isOK,:);
    end
end
    
% sets/calculates the raw, mean and SEM proportional activity values
for i = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf('%s (Apparatus %i of %i)',wStr{2},i,nApp);
    if h.Update(2,wStrNw,(2*i-1)/(2*nApp+1))
        [plotD,ok] = deal([],false);
        return
    end    
    
    % memory allocation
    IBIR = cell(nExp,nFly(i));
    
    % sets all of the experiments into a single array
    for j = 1:nExp
        if ~isempty(I{i,j})
            % sets the column indices and increments the offset       
            I{i,j}(1,:) = false;            
            
            % groups the frames where the fly has moved. from this,
            % calculate the duration of these movement events.
            moveGrp = cellfun(@(x)(getGroupIndex(x)),...
                            num2cell(I{i,j},1),'un',0);
            stopGrp = cellfun(@(x)(getGroupIndex(~x)),...
                            num2cell(I{i,j},1),'un',0);          
                        
            % removes the first group (remove stationary flies and the
            % first frame which is counted as being an inactive frame)
            stopGrp = cellfun(@(x)(x(2:end)),stopGrp,'un',0);             
            if ~isempty(indT{j})
                for k = 1:length(stopGrp)
                    % determines the stop time groups that encompass the
                    % large time gaps and removes them                    
                    indX = cell2mat(arrayfun(@(x)(find(cellfun(@(y)...
                             (any(y==x)),stopGrp{k}))),indT{j},'un',0));
                    if ~isempty(indX)
                        B = ~setGroup(indX,size(stopGrp{k}));
                        stopGrp{k} = stopGrp{k}(B);
                    end
                end
            end            
            
            % calculates the duration of the movement/stopped events
            tGrpMove = cellfun(@(x)(cellfun(@(y)(TT{j}(y(end)) - ...
                            TT{j}(y(1)-1)),x)),moveGrp,'un',0);
            tGrpStop = cellfun(@(x)(cellfun(@(y)(TT{j}(y(end)) - ...
                            TT{j}(y(1)-1)),x)),stopGrp,'un',0);                     
                        
            % determine which groups had a movement duration greater than
            % the time, tMove. these events are considered to be start
            % events. sum up the number of starts for each fly
            isStart = cellfun(@(x)(x > cP.tMove),tGrpMove,'un',0);
            
            % removes the non-starting values from the moving time groups
            % and places them in the stopped time groups
            tGrpStop = cellfun(@(x,y,z)([x;y(~z)]),...
                            tGrpStop,tGrpMove,isStart,'un',0);                       
                
            % sets the raw action initiation values
            kk = ~cellfun('isempty',tGrpStop);
            ifok = find(snTot(j).iMov.flyok{i});
            IBIR(j,ifok(kk)) = tGrpStop(kk);
        end        
    end
    
    % updates the waitbar figure
    wStrNw = sprintf(...
            'Weibull Statistics Calculations (Apparatus %i of %i)',i,nApp);
    h.Update(2,wStrNw,2*i/(2*nApp+1));
    a = calcIBIStats(IBIR);
        
    % sets the mean fitted parameters and signal data
    [plotD(i).kmn,plotD(i).lmn] = deal(a.mWP.k,a.mWP.l);
    [plotD(i).xB,plotD(i).YmnR] = deal(a.mWP.xB,a.mWP.yR);
    plotD(i).YmnF = a.mWP.yF;
        
    % sets the individual parameter values and metrics
    szR = size(plotD(i).k);
    if length(szR) == 2; szR = [szR,1]; end
    [k,l,B,M] = field2cell(a.iWP,{'k','l','B','M'},1);  
    
    % sets the shape of the arrays based on the current size
    szRnw = [szR(1), min(szR(2), size(k,2)), min(szR(3), size(k,1))];
    [jj, kk] = deal(1:szRnw(2), 1:szRnw(3));
    
    % stores the values in the plotting data struct
    plotD(i).k(:, jj, kk) = num2cell(reshape(k',szRnw));
    plotD(i).l(:, jj, kk) = num2cell(reshape(l',szRnw));
    plotD(i).B(:, jj, kk) = num2cell(reshape(B',szRnw));
    plotD(i).M(:, jj, kk) = num2cell(reshape(M',szRnw));
            
    % calculates the metric statistics
    plotD(i) = calcMetricStats(plotD(i),'k'); 
    plotD(i) = calcMetricStats(plotD(i),'l');
    plotD(i) = calcMetricStats(plotD(i),'B');    
    plotD(i) = calcMetricStats(plotD(i),'M');                
end
    
% closes the waitbar figure
if ~h.Update(2,'Activity Calculations Complete!',1)
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
[p,hLg] = deal(plotD{1}(ind),[]);

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% retrieves the formatting struct
pF = retFormatStruct(pF,1);
xi = 1:length(ind);

% sets the title string
pF.Title.String = pP.pMet;

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% sets up the plot axis
hAx = createSubPlotAxes();    
axis(hAx,'on'); hold on; 

% sets the legend strings
pF.Legend.String = pData.appName(ind);

% sets the plot values for the 
if (strcmp(pP.pMet,'Mean Survival Curve'))
    % retrieves the plot values
    [Yplt,YpltF,Xplt] = field2cell(p,{'YmnR','YmnF','xB'});         
    pF.xLabel.String = 'Time (sec)';
    pF.yLabel.String = 'Probability';   
    
    % retrieves the colours
    col = num2cell(distinguishable_colors(nApp,'w'),2);
    
    % plots the survival curves
    hPlot = cell(nApp,1);
    for i = 1:nApp
        hPlot{i} = plot(hAx,Xplt{i},Yplt{i},'color',col{i},'LineWidth',pP.pL);
        if (pP.plotFit)
            % plots the fitted values (if required)
            plot(hAx,Xplt{i},YpltF{i},'--','Color',col{i},'LineWidth',pP.pL);
        end
    end
            
    % creates the legend object
    hLg = createLegendObj(hPlot,pF.Legend,1,0);
    
    % updates the axis properties
    set(hAx,'yscale','log','box','on');   
elseif (strcmp(pP.pMet,'Burstiness + Memory'))
    % retrieves the plot values
    [Xplt,Yplt] = field2cell(p,{'M','B'});         
    pF.xLabel.String = 'Memory (M)';
    pF.yLabel.String = 'Burstiness (B)';   
    
    % retrieves the colours
    col = num2cell(distinguishable_colors(nApp,'w'),2);    
    
    % plots the survival curves
    hPlot = cell(nApp,1);
    for i = 1:nApp
        [XpltNw,YpltNw] = deal(cell2mat(Xplt{i}(:)),cell2mat(Yplt{i}(:)));
        hPlot{i} = plot(hAx,XpltNw,YpltNw,'o','color',col{i});
    end    
    
    % creates the legend object
    hLg = createLegendObj(hPlot,pF.Legend,1,0);
    
    % updates the axis properties
    set(hAx,'box','on');      
else
    % other initialisations
    xLim = xi([1 end])+0.5*[-1 1];    
    
    % sets the plot variable string and y-axis string
    switch (pP.pMet)        
        case ('Shape Factor')
            [pF.yLabel.String,pStr] = deal('Shape (k)','k');        
        case ('Scale Factor')
            [pF.yLabel.String,pStr] = deal('Shape (l)','l');
        case ('Burstiness')
            [pF.yLabel.String,pStr] = deal('Burstiness (B)','B');
        case ('Memory')
            [pF.yLabel.String,pStr] = deal('Memory (M)','M');
    end
    
    % plots the bar/boxplot graph
    plotBarBoxMetrics(hAx,xi,p,pStr,pP,[],'b');
    
    % updates the axis properties
    set(hAx,'xticklabel',[],'xlim',xLim,'linewidth',1.5,'UserData',1);
end
    
% % sets the plot
% if (strcmp(pP.pMet,'Inter-Bout Interval'))
%     lyMax = ceil(log10(max(get(hAx,'ylim'))));
%     set(hAx,'ylim',[10^(min(lyMax-1,-1)) 10^lyMax])
% end
      
% adds in the gridlines (if checked)
if (pP.plotGrid); grid(hAx,'on'); end

% ------------------------------ %
% --- PLOT AXES REFORMATTING --- %
% ------------------------------ %

% sets the x-axis labels
formatPlotAxis(hAx,pF,1);

% formats and resets the axis positions
resetAxesPos(hAx,1,1); 

% sets the group strings (if not plotting the survival curve)
if (~isempty(hLg))
    % creates the legend object       
    lgP = resetVertLegendWidth(hLg);
    lgP(1:2) = [(1-lgP(3)),(0.5-lgP(4)/2)];
    set(hLg,'position',lgP);      

    % resets the width of the trace plot axis
    pAx = get(hAx,'position');
    resetObjPos(hAx,'width',lgP(1) - (pAx(1)+0.025))    
else
    setGroupString(hAx,pF,xi,snTot(1).iMov.pInfo.gName(ind),30);
end
    
% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% ----------------------------------------------------------------------- %
% ---                         OTHER FUNCTIONS                         --- %
% ----------------------------------------------------------------------- %

% --- calculates the IBI statistics
function ibiStats = calcIBIStats(IBI)

% memory allocation and other initialisations
gof = cell(size(IBI));
[k,l,B,M] = deal(NaN(size(IBI)));
g = fittype('A+k*x','coeff',{'A','k'});

% sets the ibi statistics parameter struct
ibiStats = struct('mWP',[],'iWP',[]);

% calculates the parameters for all flies
[kT,lT,gofT,yTR,xT0] = fitWeibullPara(cell2mat(IBI(:)),g);
if (isnan(kT))
    return; 
elseif (gofT.rsquare < 0.5)
    return; 
end

% sets the mean data values into the data struct
[yF,yTR(1)] = deal(exp(-(xT0/lT).^kT),1);
mWP = struct('k',kT,'l',lT,'gof',gofT,'xB',xT0,'yR',yTR,'yF',yF);

% calculates the individual survival curve weibull parameters
y = NaN(length(xT0),numel(IBI));
for i = 1:numel(IBI)        
    % calculates the weibull parameters
    if ~isempty(IBI{i})
        N = length(IBI{i});
        [k(i),l(i),gof{i},y(:,i)] = fitWeibullPara(IBI{i},g,xT0);

        % calculates the other statistics (if feasible
        if ~isnan(k(i))
            % calculates the burstiness values 
            mIBI = mean(IBI{i},'omitnan');
            sIBI = std(IBI{i},[],'omitnan');
            B(i) = (sIBI - mIBI)/(sIBI + mIBI);

            % calculates the memory values
            [t1,t2] = deal(IBI{i}(1:(N-1)),IBI{i}(2:N));
            M(i) = sum((t1-mean(t1)).*(t2-mean(t2)))/((N-1)*std(t1)*std(t2));
        end
    end
end

% sets the individual fly statistics
iWP = struct('k',[],'l',[],'gof',[],'B',[],'M',[]);
[iWP.k,iWP.l,iWP.gof,iWP.B,iWP.M] = deal(k,l,gof,B,M);

% returns the overall data struct
[ibiStats.mWP,ibiStats.iWP] = deal(mWP,iWP);

% ---
function [k,l,gof,y,x0] = fitWeibullPara(IBI,g,x0)

%
nMin = 20;

% determines if there are any feasible data values
if (length(IBI) > nMin)
    % determines if the bin-size has already been determined
    if (nargin < 3)
        % if not, then determine the bin size from the dataset
        bSz = max(min(diff(unique(IBI))),max(IBI)/1000);
        IBI = IBI - min(IBI);
    
        % Calculate histogram
        x0 = (0:bSz:max(IBI))';
    else
        % otherwise, calculate the bin size from the previous data
        bSz = diff(x0([1 2]));
    end
       
    % calculates the IBI time histogram
    y = hist(IBI,x0);
    y = y(:)/sum(y,'omitnan');
    
    % Calculate surivival histogram
    x = x0 - bSz/2;
    [x(1),y(end:-1:1)] = deal(0,cumsum(y(end:-1:1)));
        
    % Linearize
    [y(~logical(y)),y(y>0.9999)] = deal(NaN);
    [xL,yL] = deal(log(x),log(-log(y)));
    ii = ~(isinf(xL) | isnan(xL) | isinf(yL) | isnan(yL));
    
    % Calculate linear fit
    if sum(ii) >= 2 
        % calculates the fit parameters (if enough fit points)
        fOpt = fitoptions('Method','NonlinearLeastSquares',...
                          'StartPoint',rand(1,2));
        [pExp,gof] = fit(xL(ii),yL(ii),g,fOpt);
        [k,l] = deal(pExp.k,exp(-pExp.A/pExp.k));    
    else
        % not enough feasible values, so return NaN/empty arrays
        [k,l,y,gof] = deal(NaN,NaN,NaN(size(x0)),[]);
    end
else
    % no feasible values, so return NaN/empty arrays
    [k,l,y,gof] = deal(NaN,NaN,NaN(size(x0)),[]);
end