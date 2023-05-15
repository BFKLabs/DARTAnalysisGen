% --- calculates the reaction proportion of the flies with respect to 
%     pre-stimuli sleep duration. Also calculates the fly activity 
%     exponential response metrics (long experiment only)
function pData = SleepIntensity(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Sleep Intensity (Stimuli Response)';
pData.Type = {'Pop','Multi'};
pData.fType = [2 2 3 1];
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
    [pData.hasSR,pData.hasSP,pData.hasRS] = deal(true,true,false);
    pData.sP = initSpecialPara(snTotL,pData,'nBin',1); 
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
nPara = 8;                       
cP = setParaFields(nPara);

% sets the parameter list for 
pList = cellfun(@num2str,num2cell([4:6 10 12 15 20 30]),'un',0);
pList2 = {'Mean','Median'};

% tab list headings
a = {'1 - General','2 - Stimuli Response'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'List',{4,pList},'nBin','Grouped Time Bin Size (min)');
cP(2) = setParaFields(a{1},'Boolean',1,'sepDN','Separate Results By Day/Night');
cP(3) = setParaFields(a{2},'List',{1,pList2},'cType','Speed Calculation Type');
cP(4) = setParaFields(a{2},'Number',2,'tBefore','Pre-Stimuli Signal Duration (min)',[1 10 true]);
cP(5) = setParaFields(a{2},'Number',15,'tAfter','Post-Stimuli Signal Duration (min)',[1 45 true]);
cP(6) = setParaFields(a{2},'Boolean',0,'useDouble','Use Double Exponential Inactivity Fit');
cP(7) = setParaFields(a{2},'Boolean',0,'useToff','Use Exponential Time Offset');
cP(8) = setParaFields(a{2},'Boolean',0,'moveOnly','Ignore Non-Reactive Flies For Stimuli Response');
cP(9) = setParaFields(a{2},'Boolean',0,'ignoreWeak','Ignore Weak Signal Response');
cP(10) = setParaFields(a{2},'Number',0.001,'pWeak','Peak/Signal Length Ratio', [0 0.1 true],{9,2});

% sets the tool-tip strings
cP(1).TTstr = 'Duration of the immobility time discretisation bins';
cP(2).TTstr = 'Seperates activity into day and night time groupings';
cP(3).TTstr = 'The speed calculation type (mean or median)';
cP(4).TTstr = 'The signal duration before the stimuli event';
cP(5).TTstr = 'The signal duration after the stimuli event';
cP(6).TTstr = 'Uses double exponential fit instead of single exponential';
cP(7).TTstr = 'Removes non-reactive flies from the stimuli response signal';
cP(8).TTstr = 'Duration for which a fly is considered to be non-reactive';
cP(9).TTstr = 'Ignore signals that are too weak to be analysed';
cP(10).TTstr = 'Peak/Signal length ratio beneath which the signal is considered too weak';

% adds the unique motor parameters
cP = addUniqueMotorPara(cP,snTot);

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 10;
pP = setParaFields(nPara);

% sets the general parameter lists
pList = {'Stimuli Response','Sleep Intensity Histograms'};

% sets the stimuli response parameter lists
pListGOF = {'R-Squared','Adjusted R-Squared'};

% sets the sleep intensity parameter lists
pListH = {'Histogram Count','Reaction Proportion','Combined Histograms'};

% tab list headings
a = {'1 - General','2 - Stimuli Response','3 - Histograms'};

% sets the parameter fields
pP(1) = setParaFields(a{1},'List',{1,pList},'pMet','Figure Plot Type',[],{{{'SR'},{'SP'}},{1,2}});
pP(2) = setParaFields(a{1},'Number',0.9,'pW','Metric Bar Plot Relative Width',[0 1 false]);
pP(3) = setParaFields(a{1},'Boolean',0,'plotGrid','Show Axes Gridlines');
pP(4) = setParaFields(a{2},'List',{1,pListGOF},'gofType','Goodness-of-fit Type',[],{1,1});
pP(5) = setParaFields(a{2},'Boolean',0,'plotErrorS','Plot Stimuli Response SEM',[],{1,1});
pP(6) = setParaFields(a{2},'Boolean',1,'plotRaw','Plot Raw Trace Stimuli Response Amplitude',[],{1,1});
pP(7) = setParaFields(a{2},'Boolean',1,'plotFixedS','Fix Signal Limit To Overall Maximum',[],{1,1});
pP(8) = setParaFields(a{2},'Boolean',1,'plotFixedH','Fix Metric Limit To Overall Maximum',[],{1,1});
pP(9) = setParaFields(a{2},'Boolean',1,'showGOF','Show GOF Statistics',[],{1,1});
pP(10) = setParaFields(a{3},'List',{1,pListH},'pMetH','Histogram Plot Metrics',[],{1,2});
pP(11) = setParaFields(a{3},'Boolean',0,'grpType','Group Histograms By Type',[],{1,2});

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot,Type)

% sets the figure type
if (nargin == 1); Type = 1; end

% initialises the plot format struct based on the type
if (Type == 1)
    % memory allocation
    nSub = 6;
    pF = setFormatFields(nSub);

    % initialises the font structs
    pF.Title = setFormatFields(setupFontStruct('FontSize',20),'',nSub);
    pF.xLabel = setFormatFields(setupFontStruct('FontSize',14),'',nSub);
    pF.yLabel = setFormatFields(setupFontStruct('FontSize',14),'',nSub);
    pF.Axis = setFormatFields(setupFontStruct('FontSize',10),[]);

    % sets the apparatus names as the titles
    pF.xLabel(1).String = 'Time Relative To Stimuli Event (min)';
    for i = 2:6; pF.xLabel(i).String = ''; end

    % sets the apparatus names as the titles
    pF.yLabel(1).String = 'Speed (mm s^{-1})';
    pF.yLabel(2).String = 'R^{2}';
    pF.yLabel(3).String = 'Speed (mm s^{-1})';
    pF.yLabel(4).String = 'Time (min)';
    pF.yLabel(5).String = 'Time (min)';
    pF.yLabel(6).String = 'Time (min)';

    % sets the apparatus names as the titles
    pF.Title(1).String = 'Stimuli Response';
    pF.Title(2).String = 'GOF Stats';
    pF.Title(3).String = 'Amplitude';
    pF.Title(4).String = '\tau_{Activation}';
    pF.Title(5).String = '\tau_{Inactivation}';
    pF.Title(6).String = '\tau_{Slow}';

    % sets the apparatus names as the titles
    pF.Axis(1).Font.FontSize = 14;
    pF.xLabel(1).Font.FontSize = 16;
    pF.yLabel(1).Font.FontSize = 16;
else
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
end
    
% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);
Stats1 = {'ZTestGroup','Tgrp','dnStr'};
Stats2 = {'TTestGroup','Tgrp','dnStr'};
Stats3 = {'GOF','Tgrp','dnStr'};
[xDep1,xDep2] = deal({'T','Tgrp','dnStr'},{'Tgrp','dnStr'});
[Type1,Type2,Type3] = deal(4,3,2);

% sets the independent variable fields
oP = addXVarField(oP,'Time','T','Time');
oP = addXVarField(oP,'Time Bin','Tgrp','Group');
oP = addXVarField(oP,'Day/Night','dnStr','Other');

% sets the dependent variable fields
oP = addYVarField(oP,'Speed (Abs)','Y',[],Type1,xDep1);
oP = addYVarField(oP,'Speed (Rel)','Y_rel',[],Type1,xDep1);
oP = addYVarField(oP,'Speed (SEM)','Y_sem',[],Type1,xDep1);
oP = addYVarField(oP,'Speed (Fitted)','Y_fit',[],Type1,xDep1);
oP = addYVarField(oP,'Histogram Count','Hist',[],Type2,xDep2,1);
oP = addYVarField(oP,'Reaction Count','HistR',[],Type2,xDep2,1);
oP = addYVarField(oP,'Histogram Count','Pr_N',[],Type3,xDep2);
oP = addYVarField(oP,'Reaction Proportion','Pr',Stats1,[],xDep2);
oP = addYVarField(oP,'Reaction Proportion (Mean)','Pr_P',[],Type3,xDep2);
oP = addYVarField(oP,'Reaction Proportion (SEM)','Pr_sem',[],Type3,xDep2);
oP = addYVarField(oP,'Fitted Amp','Yamp',Stats2,[],xDep2);
oP = addYVarField(oP,'Fitted Amp (Mean)','Yamp_mn',[],Type3,xDep2);
oP = addYVarField(oP,'Fitted Amp (SEM)','Yamp_sem',[],Type3,xDep2);
oP = addYVarField(oP,'Raw Amp','YampR',Stats2,[],xDep2);
oP = addYVarField(oP,'Raw Amp (Mean)','YampR_mn',[],Type3,xDep2);
oP = addYVarField(oP,'Raw Amp (SEM)','YampR_sem',[],Type3,xDep2);
oP = addYVarField(oP,'Pre-Stim Spd','YampR',Stats2,[],xDep2);
oP = addYVarField(oP,'Pre-Stim Spd (Mean)','YampR_mn',[],Type3,xDep2);
oP = addYVarField(oP,'Pre-Stim Spd (SEM)','YampR_sem',[],Type3,xDep2);
oP = addYVarField(oP,'Act Tau','kA',Stats2,[],xDep2);
oP = addYVarField(oP,'Act Tau (Mean)','kA_mn',[],Type3,xDep2);
oP = addYVarField(oP,'Act Tau (SEM)','kA_sem',[],Type3,xDep2);
oP = addYVarField(oP,'Inact Tau 1','kI1',Stats2,[],xDep2);
oP = addYVarField(oP,'Inact Tau 1 (Mean)','kI1_mn',[],Type3,xDep2);
oP = addYVarField(oP,'Inact Tau 1 (SEM)','kI1_sem',[],Type3,xDep2);                        
oP = addYVarField(oP,'Inact Tau 2','kI2',Stats2,[],xDep2);
oP = addYVarField(oP,'Inact Tau 2 (Mean)','kI2_mn',[],Type3,xDep2);
oP = addYVarField(oP,'Inact Tau 2 (SEM)','kI2_sem',[],Type3,xDep2);
oP = addYVarField(oP,'Max Response Time','Tmax',[],Type3,xDep2);
oP = addYVarField(oP,'Signal GOF','gof',Stats3,[],xDep2);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
Tgrp = dcObj.plotD{1}(1).Tgrp;
pP = retParaStruct(dcObj.pData.pP);
cP = retParaStruct(dcObj.pData.cP);
sType = {'Day Activity','Night Activity'};
[tStr,sStr] = deal('Bin Size','Temporal Activity');

% sets the specific metric class fields
switch pP.pMet
    case 'Stimuli Response'
        % case is the stimuli response metrics
        
        % field retrievals
        iAx = dcObj.getSelectAxesIndex;
        uStr = {'mm/sec','unitless','mm/sec','min','min'};
        mStr = {'Speed','GOF','Amplitude','Activation TC','Inactivation TC'};

        % removes the GOF field (if not being shown)
        if ~pP.showGOF
            ii = ~setGroup(2,size(uStr));
            [uStr,mStr] = deal(uStr(ii),mStr(ii));
        end
        
        % incorporates the fields for the double-exponential field
        if cP.useDouble
            mStr{end} = 'Fast Inactivation TC';
            [uStr{end+1},mStr{end+1}] = deal('min','Slow Inactivation TC');
        end
        
        % sets the common class fields
        dcObj.yName = mStr{iAx};
        dcObj.yUnits = uStr{iAx};
        dcObj.grpName = dcObj.pData.appName;
        
        % sets up the axes specfic fields
        if iAx == 1
            % case is the average speed trace

            % sets the histogram class fields
            dcObj.xName = 'Time';
            dcObj.xName2 = 'Bin Size';            
            dcObj.xUnits = 'sec';
            dcObj.xUnits2 = 'min';            
            dcObj.pType = 'Multi-Fitted Trace'; 
            [dcObj.xGrp,dcObj.xGrp2] = deal(Tgrp,sType);
            
        else
            % case is the other metrics

            % sets the histogram class fields
            dcObj.xGrp = sType;
            dcObj.xName = 'Time Bin';
            dcObj.pType = 'Multi-Bar Graph';

            % case is the data is grouped by time-groups
            [dcObj.xName,dcObj.xName2] = deal(sStr,tStr);
            [dcObj.xGrp,dcObj.xGrp2] = deal(sType,Tgrp);            
            
        end
        
    case 'Sleep Intensity Histograms'
        % case is the sleep intensity histograms
        
        % sets the histogram class fields
        dcObj.yName = pP.pMetH;
        dcObj.xName = 'Time Bin';
        dcObj.pType = 'Multi-Bar Graph';
        
        % sets the y-metric units
        switch pP.pMetH
            case 'Histogram Count'
                % case is the histogram count
                dcObj.yUnits = 'count';
                
            case 'Reaction Proportion'
                % case is the reaction proportion
                dcObj.yUnits = '%';
        end
        
        % initialisations
        grpName = dcObj.pData.appName;        
        
        % sets the x-metric fields (based on grouping type)
        if pP.grpType
            % case is the data is grouped by day/night activity
            dcObj.grpName = sType;
            [dcObj.xName,dcObj.xName2] = deal(tStr,'Group Name');
            [dcObj.xGrp,dcObj.xGrp2] = deal(Tgrp,grpName);
            
        else
            % case is the data is grouped by time-groups
            dcObj.grpName = grpName;
            [dcObj.xName,dcObj.xName2] = deal(tStr,sStr);
            [dcObj.xGrp,dcObj.xGrp2] = deal(Tgrp,sType);
        end
        
end

% sets up the data cursor string
dTxt = dcObj.setupCursorString();

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
tBin = str2double(cP.nBin);
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);

% memory allocation
[nGrp,tBefore] = deal(60/tBin,cP.tBefore*60);
[Tbin,cP.nGrp] = deal((tBin:tBin:60)-(tBin/2),num2str(nGrp));

% setting up of the stimuli signal time vector
T = setStimTimeVector(cP);

% retrieves the other calculation parameters (if they exist)
[dT,chT] = deal([]);
if isfield(cP,'devType'); dT = cP.devType; end
if isfield(cP,'chType'); chT = cP.chType; end

% ensures to include time bins of adequate length
Ts0 = arrayfun(@(x)(getMotorFiringTimes(x.stimP,dT,chT)),snTot,'un',0);
dTs = cellfun(@(x)(min(diff(x))),Ts0);
Tbin = Tbin(Tbin <= max(convertTime(dTs,'sec','min')));

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

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                             'T',T','Tgrp',Tgrp,'Tbin',Tbin','dnStr',dnStr,...
                             'Y',[],'Y_rel',[],'Y_sem',[],'Y_fit',[],...
                             'Yamp_mn',[],'YampR_mn',[],'Y0_mn',[],...
                             'Yamp_sem',[],'YampR_sem',[],'Y0_sem',[],...
                             'kA_mn',[],'kA_sem',[],'kI1_mn',[],'kI1_sem',[],...
                             'kI2_mn',[],'kI2_sem',[],'Tmax',[],'gof',[],...                             
                             'Pr_P',[],'Pr_N',[],'Pr_sem',[],...
                             'Hist',[],'HistR',[],'indCombMet','sum');                                  

% initialises the raw data arrays
plotD = initRawArray(plotD,NaN(cP.sepDN+1,nGrp),{'Hist','HistR'});                          
                         
% other initialisations 
is2D = ~isempty(snTot(1).Py);  

% other memory allocations
P = repmat({cell(nGrp,1+cP.sepDN)},nApp,1+is2D);
Ycount = repmat({zeros(1+cP.sepDN,nGrp)},nApp,1);

% ---------------------------------------------------- %
% --- STIMULI TRACE & IMMOBILITY TIME CALCULATIONS --- %
% ---------------------------------------------------- %

% sets the waitbar offset (is > 0 for more than one 
wStr = {'Overall Progress','Fitting Signal Exponentials'};
wOfs = length(snTot) > 1; 

% creates the waitbar figure
wStr = wStr((2-wOfs):end);
h = ProgBar(wStr,'Sleep Intensity Calculations');

% loops through each of the 
for i = 1:nExp
    % updates the waitbar figure (if more than one solution file)
    if wOfs > 0    
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
        h.Update(1,wStrNw,i/(1+nExp));
    end
    
    % calculates the new sleep intensity data
    [YcountNw,YcountRNw,XbinNw,YbinNw,ok] = ...
                            getStimuliResponseData(snTot(i),cP,h,wOfs);        
    if ~ok
        % if the user cancelled, then exit the function
        plotD = [];
        return
    else
        % otherwise, append the data to the arrays        
        for j = 1:nApp
            % appends the x-positional data to the arrays
            Nc = [1 1 numel(YcountNw{j})];
            if ((Nc(3) > 0) && (~isempty(XbinNw)))
                XbinNwT = cell2cell(cellfun(@(y)(cellfun(@(x)(cell2mat(x)),...
                                num2cell(cell2cell(y,0),2),'un',0)),...
                                num2cell(XbinNw{j},1),'un',0),0);
                for k1 = 1:nGrp
                    for k2 = 1:(1+cP.sepDN)
                        P{j,1}{k1,k2} = [P{j,1}{k1,k2},XbinNwT{k1,k2}];
                    end
                end
            end
            
            % appends the y-positional data to the arrays
            if (Nc(3) > 0) && ~isempty(YbinNw)
                YbinNwT = cell2cell(cellfun(@(y)(cellfun(@(x)(cell2mat(x)),...
                                num2cell(cell2cell(y,0),2),'un',0)),...
                                num2cell(YbinNw{j},1),'un',0),0);
                for k1 = 1:nGrp
                    for k2 = 1:(1+cP.sepDN)
                        P{j,2}{k1,k2} = [P{j,2}{k1,k2},YbinNwT{k1,k2}];
                    end
                end
            end           
                                  
            if prod(Nc) > 0
                % sets the total/reaction counts
                [iR, iC] = deal(1:size(YcountNw{j},1), 1:size(YcountNw{j},2));
                plotD(j).Hist(iR,iC,i) = YcountNw{j};
                plotD(j).HistR(iR,iC,i) = YcountRNw{j}; 

                % adds on the bin/reaction counts                
                Ycount{j} = Ycount{j} + sum(cell2mat(reshape(YcountNw{j},Nc)),3);
%                 YcountR{j} = YcountR{j} + sum(cell2mat(reshape(YcountRNw{j},Nc)),3);
            end
        end
    end    
end

% loops through each of the apparatus calculating the exponentials
xiT = max(1,tBefore-60):tBefore;
[jj,a,N] = deal((tBefore+1):length(T),min(1,wOfs),[cP.tBefore,cP.tAfter]*60);
for i = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf('Fitting Exponential Curves (Region %i of %i)',i,nApp);
    if h.Update(1+wOfs,wStrNw,0.5*(1-a)+0.5*(1+a)*(i+1)/(2+nApp))
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    end                
    
    % calculates the metrics/signals for all groups
    for j = 1:(1+cP.sepDN)
        if ~all((Ycount{i}(j,:)) == 0)
            % calculates the average speed        
            [plotD(i).Y,plotD(i).Y_sem] = ...
                                calcSRAvgSpeed(P(i,:),cP,N,[],[],is2D,1);
                            
            % calculates the initial speed and the relative/absolute 
            % stimuli response speed values   
            plotD(i).Y0_mn = cellfun(@(x)(median(x(xiT,:),1,'omitnan')),...
                                plotD(i).Y,'un',0);        
            plotD(i).Y_rel = cellfun(@(x,y)(x-repmat(y,size(x,1),1)),...
                                plotD(i).Y,plotD(i).Y0_mn,'un',0);
            
            % fits exponentials to each of the time-binned groups   
            Ydata = cellfun(@(x)(x(jj,:)),plotD(i).Y_rel,'un',0);   
            [p,Yfit,gof] = fitSignalExp(plotD(i).T(jj),Ydata,cP);
            plotD(i).gof = cellfun(@(x)(x(:)'),gof,'un',0);
            plotD(i).Y_fit = cellfun(@(x)....
                                ([zeros(tBefore,nGrp);x]),Yfit,'un',0);            
                            
            % sets the SR fitted parameters
            plotD(i) = setSRFittedPara(plotD(i),p,tBefore,nGrp);                                    
        end
    end
    
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
    
%     % sets the histogram counts and reaction proportion
%     [plotD(i).Pr_N,noCount] = deal(Ycount{i}, Ycount{i}==0);
%     [plotD(i).Pr_P,PR] = deal(YcountR{i}./Ycount{i});
%     plotD(i).Pr_sem = sqrt((PR.*(1 - PR))./Ycount{i});
%     [plotD(i).Pr_P(noCount),plotD(i).Pr_sem(noCount)] = deal(0);
end
   
% closes the waitbar figure
if ~h.Update(1,'Inactivity Calculations Complete!',1)
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

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the title string
grpStr = snTot(1).iMov.pInfo.gName{sP.pInd};
pF.Title(1).String = sprintf('%s (%s)',pF.Title(1).String,grpStr);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% creates the figure based on the type
switch pP.pMet
    case ('Stimuli Response')
        pF = initPlotFormat([],1);
        createStimResFig(hP,plotD{1},cP,pP,sP,pF);
    case ('Sleep Intensity Histograms')
        pF = initPlotFormat(snTot(1),2);
        createSIHistFig(hP,pData,plotD{1},cP,pP,sP,pF);
end
  
% --- creates the sleep intensity histogram figures
function createSIHistFig(hP,pData,plotD,cP,pP,sP,pF)

% determines if calculates are separated by day/night
sepDN = size(plotD(1).Pr_P,1) == 2;

% if the incorrect combination is used, then exit with an error (not
% possible to plot the day/night separation and combined histograms
% together)
if strcmp(pP.pMetH,'Combined Histograms') && sepDN
    eStr = 'Not possible to plot the Combined Histograms with Day/Night Separation';
    waitfor(msgbox(eStr,'Incorrect Plot Format','modal'))
    return
elseif strcmp(pP.pMetH,'Combined Histograms') && pP.grpType
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
p = plotD(ind);

% sets up the x-axis strings
xStr = p(1).Tgrp;

% other parameters
col = 'yk';

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% retrieves the formatting struct
if isempty(m)
    szMx = 1;
elseif pP.grpType
    szMx = ceil(length(sP.Sub.isPlot)/2);
    pF.Legend.String = pData.appName(ind);
    [ind,n] = deal(1:(1+cP.sepDN),1);
else
    szMx = max([m n]);
end

% sets the y-label string
pF = retFormatStruct(pF,szMx);

% sets the x/y-label strings
pF.xLabel.String = 'Immobility Time (min)'; 
pF.yLabel.String = 'Frequency'; 

% removes the label indices
[pF.xLabel.ind,pF.yLabel.ind,pF.zLabel.ind] = deal(NaN);

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% plots the histograms for each of the apparatus
if pP.grpType
    % memory allocation
    hAx = cell(1+sepDN,1);
    
    % sets the title strings based on the calculation type
    if sepDN
        % case is day/night separation
        pF.Title = repmat(pF.Title,[1,2]);
        pF.Title(1).String = 'Day Activity';
        pF.Title(2).String = 'Night Activity';
    else
        % case is day/night combined
        pF.Title(1).String = 'Day/Night Combined Activity';
    end
        
    % sets the histogram arrays            
    N = field2cell(p,'Pr_N');        
    
    %    
    for i = 1:(1+sepDN)
        % memory allocation
        hAx{i} = createSubPlotAxes(hP,[(1+sepDN),1],i); 
        hold(hAx{i},'on');
        HistNw = cell2mat(cellfun(@(x)(x(i,:)),N,'un',0))';                

        switch pP.pMetH
            case ('Histogram Count') 
                [yLim,Yplt,Ysem] = deal([0 detOverallLimit(N)],HistNw,[]);   
            case ('Reaction Proportion')                    
                Yplt = cell2mat(cellfun(@(x)(x(i,:)),...
                                field2cell(p,'Pr_P'),'un',0))';
                [yLim,Ysem] = deal([0 100],100*sqrt((Yplt.*(1 - Yplt))./HistNw));
                Yplt = Yplt*100;
                pF.yLabel.String = '% Active'; 
        end                        
        
        % creates the bar plot
        hBar = plotBarError(hAx{i},Yplt,Ysem,~isempty(Ysem));             
        set(hAx{i},'ylim',yLim,'box','on')
        
        % formats the plot axis
        formatPlotAxis(hAx{i},pF,i,1);       
         
        % turns the grid on (if specified)
        if pP.plotGrid; set(hAx{i},'ygrid','on','yminorgrid','on'); end                       
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
        setGroupString(hAx{i}(1),pF,p(1).Tbin,xStr,-45);
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
        switch (pP.pMetH)
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
        if ~isempty(Ysem)              
            for j = 1:size(Yplt,1)
                Yplt(isnan(Ysem)|(Ysem==0)) = NaN;
                xiNw = xi + sepDN*(2*(j-1)-1)/4;
                addBarError(hAx{i}(1),xiNw,Yplt(j,:),Ysem(j,:),'r',3);
            end
        end    

        % updates the axis properties                
        if length(hAx{i}) == 1
            % formats the single axis
            set(gcf,'CurrentAxes',hAx{i}(1))
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
        if (strcmp(pP.pMetH,'Combined Histograms'))
            set(hAx{i}(2),'Position',get(hAx{i}(1),'position'))
        end
    end     
end

% --- creates the stimuli response figure
function createStimResFig(hP,plotD,cP,pP,sP,pF)

% sets the indices to be plotted
[iPlotT,iPlotF] = deal(find(sP.pT),find(sP.pF));
if isempty(iPlotT) && isempty(iPlotF); return; end

% sets the number of groups and the parameter struct
p = plotD(sP.pInd);

% sets up the x-axis strings
[pY,lStyle,nGrp] = deal(0.1,{'-',':'},60/str2double(cP.nBin));
useDouble = any(cellfun(@(x)(~all(isnan(x))),p.kI2_mn));
nSub = 3 + (pP.showGOF + useDouble);

% resets the inactivation title string (if using double exponential)
if useDouble
    pF.Title(5).String = '\tau_{Fast}';
end

% sets the subplot x-tick label strings
if (cP.sepDN)
    xTickLbl = {'Day','Night'};
else
    xTickLbl = [];
end

% sets focus to the main axis
hAx = createSubPlotAxes(hP,[2,1],1);
hold(hAx,'on')
set(hAx,'OuterPosition',[0,0.5,1,1/2],'UserData',1,'box','on')

% ----------------------------------------------- %
% --- FITTED EXPONENTIAL PARAMETER BAR GRAPHS --- %
% ----------------------------------------------- %

% sets the axis maximum
pD = plotD(~cellfun('isempty',field2cell(plotD,'Y')));

% sets the parameter strings
col = getBarColourScheme(nGrp,'m');
pStr = {'gof','Yamp_mn','kA','kI1','kI2'};
if pP.plotRaw; pStr{2} = 'YampR_mn'; end
pStr = pStr((~pP.showGOF+1):(~pP.showGOF+nSub));   

% memory allocation
hAxS = cell(nSub,1);
xLim = [1 (1+cP.sepDN)]+[-0.55 0.5];

% plots the reaction proportion parameters
for i = 1:nSub    
    % creates the bar plot
    hAxS{i} = createSubPlotAxes(hP,[2,nSub],i+nSub);
    hold(hAxS{i},'on');
    set(hAxS{i},'xlim',xLim,'linewidth',1.5,'box','on','UserData',i+1) 
    
    % retrieves the stimuli response values
    [Y,Ysem,Ymx] = getSRValues(p,pD,pP,pStr{i}); 
        
    % plots the data
    [hPlot,xTick] = plotBarError(hAxS{i},Y,Ysem,~isempty(Ysem),pP.pW,col);      

    % formats the subplot axis and sets the group strings
    formatPlotAxis(hAxS{i}(1),pF,i+(1+(~pP.showGOF)));    
    if pP.showGOF && (i == 1)
        set(hAxS{i},'ylim',[0 1]);    
    else        
        if pP.plotFixedH
            resetYAxisScale(hAxS{i},Y(:),Ymx);
        else            
            resetYAxisScale(hAxS{i},Y(:));
        end
    end
    
    % sets the x-label group strings    
    xLim = [1 length(xTick)]+[-0.55 0.5];
    set(hAxS{i},'xtick',xTick,'xticklabel',xTickLbl)      
    plot(hAxS{i},[0 nGrp],[0 0],'k','linewidth',1.5)      
        
    % turns the grid on (if specified)
    if pP.plotGrid; set(hAxS{i},'ygrid','on','yminorgrid','on'); end
end

% resets the axes positions
resetAxesPos(hAxS,1,nSub,[0.01 0])

% ------------------------------- %
% --- AVERAGE VELOCITY TRACES --- %
% ------------------------------- %

% % makes the stimuli response trace axes current
% set(gcf,'CurrentAxes',hAx)

% sets the overall maximum
YSRmn = min(cell2mat(cellfun(@(x,y)(min(x{1}(:)-y{1}(:))),...
            field2cell(pD,'Y'),field2cell(pD,'Y_sem'),'un',0)));
YSRmx = max(cell2mat(cellfun(@(x,y)(max(x{1}(:)+y{1}(:))),...
            field2cell(pD,'Y'),field2cell(pD,'Y_sem'),'un',0)));

% determines the traces/fitted exponentials that are to be plotted
if cP.sepDN
    hPlotS = cell(2,1); 
else
    hPlotS = [];
end

% plots the signal traces
for i = 1:length(iPlotT)
    j = iPlotT(i);
    for k = 1:(1+cP.sepDN)
        [Tplt,Yplt] = deal(p.T/60,p.Y{k}(:,j));
        
        if pP.plotErrorS                  
            plotSignalSEM(Yplt,p.Y_sem{k}(:,j),Tplt,col{j},0.2*k)
        end        
        
        tagStr = sprintf('hRaw%i',k);
        if (i == 1) && cP.sepDN
            hPlotS{k} = plot(hAx,Tplt,Yplt,'color',col{j},...
                'linestyle',lStyle{k},'UserData',i,'Tag',tagStr);
        else
            plot(hAx,Tplt,Yplt,'color',col{j},'linestyle',lStyle{k},...
                'UserData',i,'Tag',tagStr);            
        end        
    end
end

% plots the exponential traces
for i = 1:length(iPlotF)
    j = iPlotF(i);
    for k = 1:(1+cP.sepDN)
        tagStr = sprintf('hFit%i',k);
        [Tplt,Yplt] = deal(p.T/60,p.Y_fit{k}(:,j)+p.Y0_mn{k}(j));        
        plot(hAx,Tplt,Yplt,'color',col{j},'linewidth',2,...
            'linestyle',lStyle{k},'UserData',i,'Tag',tagStr);            
    end
end

% formats the plot axis
if any(sP.pT) || any(sP.pF)
    YY = [p.Y{1}(:,iPlotT) p.Y_fit{1}(:,iPlotF)];       
    if cP.sepDN
        YY = [YY p.Y{2}(:,iPlotT) p.Y_fit{2}(:,iPlotF)];         
    end    
        
    % sets the min/max values
    if pP.plotFixedS
        [yMin,yMax] = deal(YSRmn,YSRmx);
    else
        [yMin,yMax] = deal(min(YY(:)),max(YY(:)));
    end

    % updates the y-axis limits
    delY = (yMax-yMin);        
    set(hAx,'ylim',[roundP(yMin-pY*delY,0.1) roundP(yMax+pY*delY,0.1)])    
end

% formats the plot axis
formatPlotAxis(hAx,pF,1);
if all(isnan(YY(:))); cla(hAx); return; end

% sets the axis properties (based on the plot type)
[xLim,yLim] = deal(roundP(p.T([1 end])/60),get(hAx,'ylim'));
set(hAx,'xlim',xLim,'linewidth',1.5,'UserData',1); 
plot(hAx,[0 0],yLim,'r--','linewidth',1.5,'HitTest','off')
plot(hAx,xLim,[0 0],'r--','linewidth',1.5,'HitTest','off')

% plots the gridlines (if required)
if pP.plotGrid; grid(hAx,'on'); end

% sets the new left location of the main trace plot
[pLbl,pAx] = deal(get(get(hAx,'yLabel'),'Extent'),get(hAx,'position'));
pAx(1) = pAx(3)*(xLim(1)-pLbl(1))/diff(xLim);
            
% creates the new legend    
pF.Legend.String = p(1).Tgrp;
hLg = createLegendObj(hPlot,pF.Legend,1);
resetLegendAxisPos(hAxS,hLg,[2 nSub])
resetObjPos(hLg,'bottom',0.25-retObjDimPos(hLg,4,1)/2);

% creates the new legend    
if (~isempty(hPlotS))
    pF.Legend.String = {'Day';'Night'};
    hLg2 = createLegendObj(hPlotS,pF.Legend,1);
    resetLegendAxisPos(hAx,hLg2,[2 1])
    resetObjPos(hLg2,'bottom',0.75-retObjDimPos(hLg2,4,1)/2);

    % resets the width of the trace plot axis
    lgPNw = get(hLg,'position');
    pAx(1) = retObjDimPos(hAxS{1},1,1);
    pAx(3) = lgPNw(1) - pAx(1);   
    set(hAx,'position',pAx)
else
    pAx(1) = retObjDimPos(hAxS{1},1,1);
    pAx(3) = (retObjDimPos(hLg,1,1)+retObjDimPos(hLg,3,1)) - pAx(1);
    set(hAx,'position',pAx)
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
        switch (pData.oP.yVar(jj(i)).Stats{1})
            case ('ZTestGroup')
                pData.oP.yVar(jj(i)).Stats{1} = 'ZTest';    
            case ('TTestGroup')
                pData.oP.yVar(jj(i)).Stats{1} = 'TTest';
        end
    end    
end

% sets the separation flag
pData.oP.sepGrp = cP.sepDN;