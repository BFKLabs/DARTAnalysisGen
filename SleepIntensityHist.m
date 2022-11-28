% --- calculates the reaction proportion of the flies with respect to 
%     pre-stimuli sleep duration (long experiment only)
function pData = SleepIntensityHist(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Sleep Intensity (Histograms)';
pData.Type = {'Pop','Multi'};
pData.fType = [2 2 3 1];
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
    [pData.hasSP,pData.hasRS] = deal(true,false);
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
rI.Shape = 'None';
rI.Stim = 'Motor';
rI.Spec = 'None';

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 2;                       
cP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the parameter list for the drop-down list
pList = cellfun(@num2str,num2cell([4:6 10 12 15 20 30]),'un',0);

% sets the parameter fields
cP(1) = setParaFields(a{1},'List',{4,pList},'nBin','Grouped Time Bin Size (min)');
cP(2) = setParaFields(a{1},'Boolean',1,'sepDN','Separate Results By Day/Night');

% sets the tool-tip strings
cP(1).TTstr = 'Duration of the immobility time discretisation bins';

% adds the unique motor parameters
cP = addUniqueMotorPara(cP,snTot);

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 4;
pP = setParaFields(nPara);

% sets the parameter list for 
pList = {'Histogram Count','Reaction Proportion','Combined Histograms'};

% sets the tab list names
a = {'1 - General','2 - Histograms'};

% sets the parameter fields
pP(1) = setParaFields(a{1},'Boolean',0,'plotGrid','Plot Axis Gridlines');
pP(2) = setParaFields(a{2},'List',{1,pList},'pMet','Histogram Plot Metrics');
pP(3) = setParaFields(a{2},'Boolean',0,'groupType','Group Histograms By Type');
pP(4) = setParaFields(a{2},'Number',0.75,'pW','Bar Plot Relative Width',[0 1 false]);

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
nApp = length(snTot.iMov.ok);    
pF = setFormatFields(nApp);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],'',1);
pF.yLabel = setFormatFields([],'',1);
pF.zLabel = setFormatFields([],'% Active',1);
pF.Axis = setFormatFields([],[]);

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Title(i).String = snTot.iMov.pInfo.gName{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);
[Stats,xDep] = deal({'ZTestGroup','Tgrp','dnStr'},{'Tgrp','dnStr'});
[Type1,Type2] = deal(3,2);

% sets the independent variable fields
oP = addXVarField(oP,'Time Group','Tgrp','Group');
oP = addXVarField(oP,'Stimuli Index','dnStr','Other');

% sets the dependent variable fields
oP = addYVarField(oP,'Histogram Count','Hist',[],Type1,xDep,1);
oP = addYVarField(oP,'Reaction Count','HistR',[],Type1,xDep,1);
oP = addYVarField(oP,'Histogram Count','Pr_N',[],Type2,xDep);
oP = addYVarField(oP,'Reaction Proportion','Pr',Stats,[],xDep);
oP = addYVarField(oP,'Reaction Proportion (Mean)','Pr_P',[],Type2,xDep);
oP = addYVarField(oP,'Reaction Proportion (SEM)','Pr_sem',[],Type2,xDep);

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

% memory allocation
tBin = str2double(cP.nBin);
[nGrp,nExp] = deal(60/tBin,length(snTot));

% other initialisations
ok = false;
Tbin = (tBin:tBin:60)-(tBin/2);

% sets the group strings  
Tgrp = setTimeBinStrings(tBin,nGrp);

% sets the day/night strings (if separating activity by day/night)
if cP.sepDN
    % data is separated by day/night
    dnStr = {'Day','Night'};        
else
    % no data separation
    dnStr = [];
end

% retrieves the other calculation parameters (if they exist)
if isfield(cP,'devType'); devType = cP.devType; end
if isfield(cP,'chType'); chType = cP.chType; end

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                                 'Tbin',Tbin','Tgrp',Tgrp,'dnStr',dnStr,...
                                 'Pr_P',[],'Pr_N',[],'Pr_sem',[],...
                                 'Hist',[],'HistR',[],'indCombMet','sum');

% other initialisations                             
nApp = length(snTot(1).iMov.ok);
% [Ycount,YcountR] = deal(repmat({zeros(1+cP.sepDN,nGrp)},nApp,1));
                        
% ---------------------------------------------------- %
% --- STIMULI TRACE & IMMOBILITY TIME CALCULATIONS --- %
% ---------------------------------------------------- %

% sets the waitbar offset (is > 0 for more than one 
wStr = {'Overall Progress','Fitting Signal Exponentials'};
wOfs = (length(snTot) > 1); 

% creates the waitbar figure
wStr = wStr((2-wOfs):end);
h = ProgBar(wStr,'Sleep Intensity Calculations');

% loops through each of the experiments calculating the sleep intensity
% data values
for i = 1:nExp
    % updates the waitbar figure (if more than one solution file)
    if wOfs > 0    
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
        h.Update(1,wStrNw,i/(1+nExp));
    end
    
    % calculates the new sleep intensity data
    [YcountNw,YcountRNw] = getStimuliResponseData(snTot(i),cP,h,wOfs);            
    if isempty(YcountNw)
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    else
        % otherwise, append the data to the arrays        
        for j = 1:nApp
            % sets the total/reaction counts            
            Nc = [1 1 numel(YcountNw{j})];
            if prod(Nc) > 0               
                % 
                [iR, iC] = deal(1:size(YcountNw{j},1), 1:size(YcountNw{j},2));
                plotD(j).Hist(iR,iC,i) = YcountNw{j};
                plotD(j).HistR(iR,iC,i) = YcountRNw{j};
                
%                 %
%                 ii = any(~cellfun('isempty',plotD(j).Hist(:,:,i)),2);
%                 plotD(j).Hist(:,:,i) = plotD(j).Hist(ii,:,i);
%                 plotD(j).HistR(:,:,i) = plotD(j).HistR(ii,:,i);                

%                 % adds on the bin/reaction counts
%                 Ycount{j} = Ycount{j} + sum(cell2mat(reshape(YcountNw{j},Nc)),3);
%                 YcountR{j} = YcountR{j} + sum(cell2mat(reshape(YcountRNw{j},Nc)),3);                
            end
        end
    end    
end

% sets the histogram count/reaction proportion values
for i = 1:nApp
    % calculates the reaction proportion for each fly
    R = cellfun(@(x,y)(x./y),plotD(i).HistR,plotD(i).Hist,'un',0);
    Y = cell2mat(reshape(R(:),[1 1 length(R(:))]));
    nnN = sum(~isnan(Y),3);
    
    % sets the overall histogram counts (i.e., which flies actually reacted
    % for that given time bin)
    [plotD(i).Pr_N, noCount] = deal(nnN, nnN == 0);
        
    % calculates the mean/sem reaction proportion
    [plotD(i).Pr_P,PR] = deal(mean(Y,3,'omitnan'));
    plotD(i).Pr_sem = sqrt((PR.*(1 - PR))./nnN);
    [plotD(i).Pr_P(noCount),plotD(i).Pr_sem(noCount)] = deal(0);
end

% updates the waitbar figure
if ~h.Update(1,'Sleep Intensity Data Retrieval Complete!',1)
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

% determines if calculates are separated by day/night
sepDN = size(plotD{1}(1).Pr_P,1) == 2;

% if the incorrect combination is used, then exit with an error (not
% possible to plot the day/night separation and combined histograms
% together)
if (strcmp(pP.pMet,'Combined Histograms') && (sepDN))
    eStr = 'Not possible to plot the Combined Histograms with Day/Night Separation';
    waitfor(msgbox(eStr,'Incorrect Plot Format','modal'))
    return
elseif (strcmp(pP.pMet,'Combined Histograms') && (pP.groupType))
    eStr = 'Not possible to plot the Combined Histograms with Type Grouping';
    waitfor(msgbox(eStr,'Incorrect Plot Format','modal'))
    return    
end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the plotting indices and subplot indices
[ind,m,n] = deal(find(sP.Sub.isPlot),sP.Sub.nRow,sP.Sub.nCol);
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% other parameters
col = 'yk';

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% retrieves the formatting struct
if (isempty(m))
    szMx = 1;
elseif (pP.groupType)
    szMx = ceil(length(sP.Sub.isPlot)/2);
    pF.Legend.String = pData.appName(ind);
    [ind,n] = deal(1:(1+cP.sepDN),1);
else
    szMx = max([m n]);
end

% sets the y-label string
pF = retFormatStruct(pF,szMx);

% sets the x/y-label strings
xStr = p(1).Tgrp;
pF.xLabel.String = 'Immobility Time (min)'; 
pF.yLabel.String = 'Frequency'; 

% removes the label indices
[pF.xLabel.ind,pF.yLabel.ind,pF.zLabel.ind] = deal(NaN);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% plots the histograms for each of the apparatus
if (pP.groupType)
    % memory allocation
    hAx = cell(1+sepDN,1);
    
    % sets the title strings based on the calculation type
    if (sepDN)
        % case is day/night separation
        pF.Title = repmat(pF.Title,[1,2]);
        pF.Title(1).String = 'Day Activity';
        pF.Title(2).String = 'Night Activity';
    else
        % case is day/night combined
        pF.Title(1).String = 'Day/Night Combined Activity';
    end
        
    % sets the histogram arrays            
    Hist = field2cell(p,'Pr_N');    
        
    %    
    for i = 1:(1+sepDN)
        % memory allocation
        hAx{i} = createSubPlotAxes(hP,[(1+sepDN),1],i); 
        hold(hAx{i},'on');
        HistNw = cell2mat(cellfun(@(x)(x(i,:)),Hist,'un',0))';                

        switch (pP.pMet)
            case ('Histogram Count') 
                [yLim,Yplt,Ysem] = deal([0 detOverallLimit(Hist)],HistNw,[]);   
            case ('Reaction Proportion')                    
                Yplt = cell2mat(cellfun(@(x)(x(i,:)),...
                                        field2cell(p,'Pr_P'),'un',0))';
                [yLim,Ysem] = deal([0 100],100*sqrt((Yplt.*(1 - Yplt))./HistNw));
                Yplt = Yplt*100;
                pF.yLabel.String = '% Active'; 
        end                        
        
        % creates the bar plot
        [hBar,xTick] = plotBarError(hAx{i},Yplt,Ysem,~isempty(Ysem));             
        set(hAx{i},'ylim',yLim,'box','on')
        
        % formats the plot axis
        formatPlotAxis(hAx{i},pF,i,1);       
         
        % turns the grid on (if specified)
        if (pP.plotGrid); set(hAx{i},'ygrid','on','yminorgrid','on'); end                       
    end
    
    % sets the non-aligned x/y labels
    formatMultiXYLabels(hAx,pF,[1+sepDN,1]);
    
    % formats and resets the axis positions
    resetAxesPos(hAx,(1+sepDN),1);     
    
    % creates the legend object
    if (nApp > 1)
        hLg = createLegendObj(hBar,pF.Legend);

        % resets the legend position
        lgP = get(hLg,'position');                        
        set(hLg,'position',[(1-lgP(3)),(0.5-lgP(4)/2),lgP(3:4)]);
    end

    % reformats the axis position and x-tick labels
    for i = 1:(1+sepDN)
        % sets the group strings
        setGroupString(hAx{i}(1),pF,xTick,xStr,-45);
        axis(hAx{i},'on')        
        
        % updates the axis position
        if (nApp > 1)
            axP = get(hAx{i},'position');
            axP(3) = lgP(1) - axP(1);
            set(hAx{i},'position',axP);               
        end
    end  
else
    % memory allocation and initialisations
    hAx = cell(nApp,1);
    [pF.xLabel.ind,pF.yLabel.ind,pF.zLabel.ind] = deal(NaN);
    xi = 1:length(p(1).Tbin);
        
    %
    for i = 1:nApp
        % creates a new subplot
        [Ysem,Yplt] = deal([],[]);
        hAx{i} = createSubPlotAxes(hP,[m,n],i);
        hold(hAx{i},'on');

        % creates the axis plot
        switch (pP.pMet)
            case ('Histogram Count')            
                % plots the data
                for j = 1:(1+sepDN)
                    xiNw = xi + sepDN*(2*(j-1)-1)/4;
                    bar(xiNw,p(i).Pr_N(j,:),pP.pW/(1+sepDN),col(j),'tag','hBar');    
                end
                set(hAx{i},'ylim',[0 detOverallLimit(field2cell(p,'Pr_N'))])            
            case ('Reaction Proportion')
                for j = 1:(1+sepDN)
                    xiNw = xi + sepDN*(2*(j-1)-1)/4;
                    bar(xiNw,100*p(i).Pr_P(j,:),(pP.pW)/(1+sepDN),col(j),'tag','hBar');                    

                    prY = p(i).Pr_P(j,:);
                    Ysem = [Ysem;100*sqrt((prY.*(1 - prY))./p(i).Pr_N(j,:))];                               
                    Yplt = [Yplt;100*p(i).Pr_P(j,:)];
                    pF.yLabel.String = '% Active'; 
                end                                               
                set(hAx{i},'ylim',[0 100]); 

            case ('Combined Histograms')
                % plots the double axis bar
                [Yplt,Ysem] = deal([p(i).Pr_N;100*p(i).Pr_P],[]);
                [hAx{i},hBar] = plotDoubleAxisBar(hAx{i},xi,Yplt(1,:),...
                                      Yplt(2,:),pP.pW,[],5,[NaN 100],[],...
                                      100*p(i).Pr_sem);                  
        end

        % if the SEM signal is set, then add the error bars
        if (~isempty(Ysem))                
            for j = 1:size(Yplt,1)
                Yplt(isnan(Ysem)|(Ysem==0)) = NaN;
                xiNw = xi + sepDN*(2*(j-1)-1)/4;
                addBarError(hAx{i}(1),xiNw,Yplt(j,:),Ysem(j,:),'r',3);
            end
        end    

        % updates the axis properties                
        if (length(hAx{i}) == 1)
            % formats the single axis
            formatPlotAxis(hAx{i}(1),pF,ind(i));                 
        else
            % formats the double plot axis
            pF = formatDoubleAxis(hAx{i},hBar,pF,ind(i));            
        end        
                
        % sets the x-axis strings
        if (length(hAx{i}) > 1); resetAxisPos(hAx{i}); end
                
        % turns the grid on (if specified)
        if (pP.plotGrid)
            set(hAx{i}(1),'ygrid','on','yminorgrid','on')
        end                    

        % resets the axis limits (for all axis)
        for k = 1:length(hAx{i})
            set(hAx{i}(k),'xlim',xi([1 end])+0.5*[-1 1],'box','on');    
        end        
    end
    
    % sets the non-aligned x/y labels
    formatMultiXYLabels(hAx,pF,[m,n]);    
       
    % resets the axis positions 
    resetAxesPos(hAx,m,n,[0.02 0]);     
        
    % adds the x-tick labels (and resets any double axis)
    for i = 1:nApp
        setGroupString(hAx{i}(1),pF(1),xi,xStr,-45);
        if (strcmp(pP.pMet,'Combined Histograms'))
            set(hAx{i}(2),'Position',get(hAx{i}(1),'position'))
        end
    end     
end
    
% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% initialisations
cP = retParaStruct(pData.cP);

% removes the day index dependency if not separating by day
if (~cP.sepDN)
    % determines which output variables have a day index dependency
    xDep = field2cell(pData.oP.yVar,'xDep');
    ii = find(cellfun(@(x)(any(strcmp(x,'dnStr'))),xDep));
    
    % removes the day index dependency
    for i = 1:length(ii)
        pData.oP.yVar(ii(i)).xDep = xDep{ii(i)}(~strcmp(xDep{ii(i)},'dnStr'));
    end
    
    % determines which output variables have a day index dependency
    Stats = field2cell(pData.oP.yVar,'Stats');
    jj = find(cellfun(@(x)(any(strcmp(x,'dnStr'))),Stats));
    
    % removes the day index dependency
    for i = 1:length(jj)
        pData.oP.yVar(jj(i)).Stats = Stats{jj(i)}(~strcmp(Stats{jj(i)},'dnStr'));
        pData.oP.yVar(jj(i)).Stats{1} = 'ZTest';
    end    
end

% sets the separation flag
pData.oP.sepGrp = cP.sepDN;