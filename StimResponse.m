% --- calculation of the fly activity exponential response with respect to 
%     stimuli delivered throughout the day (long experiment only)
function pData = StimResponse(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Stimuli Response (Long)';
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
    [pData.hasSR,pData.hasRS] = deal(true,false);
    pData.sP = initSpecialPara(snTotL,pData,'nGrp',1); 
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

% determines if there are any multi-day experiments
multiDay = any(detExptDayDuration(snTot) > 1);

% initialises the parameter struct
nPara = 10 + multiDay;                       
cP = setParaFields(nPara);

% sets the parameter list for the drop-down list
pList = cellfun(@num2str,num2cell([1 2 4 6 8 12]),'un',0);
pList2 = {'Mean','Median'};

% tab list headings
a = {'1 - General','2 - Exponential Fit'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'List',{2,pList},'nGrp','Number of Daily Time Groups');
cP(2) = setParaFields(a{1},'List',{1,pList2},'cType','Speed Calculation Type');
cP(3) = setParaFields(a{1},'Number',2,'tBefore','Time Before Stimuli (min)',[1 15 true]);
cP(4) = setParaFields(a{1},'Number',15,'tAfter','Time After Stimuli (min)',[1 45 true]);
cP(5) = setParaFields(a{2},'Boolean',0,'useDouble','Use Double Exponential Inactivity Fit');
cP(6) = setParaFields(a{2},'Boolean',0,'useToff','Use Exponential Time Offset');
cP(7) = setParaFields(a{1},'Boolean',1,'moveOnly','Ignore Non-Reactive Flies For Stimuli Response');
cP(8) = setParaFields(a{1},'Number',5,'tMove','Non-Reactivity Duration (min)',[1 60 true],{7,2});
cP(9) = setParaFields(a{2},'Boolean',0,'ignoreWeak','Ignore Weak Signal Response');
cP(10) = setParaFields(a{2},'Number',0.001,'pWeak','Peak/Signal Length Ratio', [0 0.1 true],{9,2});

% sets the tool-tip strings
cP(1).TTstr = 'The number of groups that the day is split up into';
cP(2).TTstr = 'The speed calculation type (mean or median)';
cP(3).TTstr = 'The signal duration before the stimuli event';
cP(4).TTstr = 'The signal duration after the stimuli event';
cP(5).TTstr = 'Uses double exponential fit instead of single exponential';
cP(6).TTstr = 'Includes a time offset for the exponential fit';
cP(7).TTstr = 'Removes non-reactive flies from the stimuli response signal';
cP(8).TTstr = 'Duration for which a fly is considered to be non-reactive';
cP(9).TTstr = 'Ignore signals that are too weak to be analysed';
cP(10).TTstr = 'Peak/Signal length ratio beneath which the signal is considered too weak';

% sets the parameter fields/tt-strings if data can be separated daily
if (multiDay)
    % sets the parameter fields
    cP(11) = setParaFields(a{1},'Boolean',0,'grpDay','Separate Response By Daily Activity');
        
    % sets the tool-tip strings
    cP(11).TTstr = 'Flag indicating whether the response is to be separated by daily activity';        
end

% adds the unique motor parameters
cP = addUniqueMotorPara(cP,snTot);

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% determines if there are any multi-day experiments
nDayMx = max(detExptDayDuration(snTot));

% initialises the parameter struct
nPara = 9 + 2*(nDayMx > 1);                      
pP = setParaFields(nPara);

% sets the parameter list for the drop-down list
pList = {'R-Squared','Adjusted R-Squared'};

% tab list headings
a = {'1 - Stimuli Trace','2 - Fit Metrics','3 - Daily Groupings'};

% sets the parameter fields
pP(1) = setParaFields(a{2},'Boolean',1,'incMet','Include Stimuli Response Metric Subplots');
pP(2) = setParaFields(a{2},'List',{1,pList},'gofType','Goodness-of-fit Type',[],{1,2});
pP(3) = setParaFields(a{1},'Boolean',1,'plotRaw','Plot Raw Trace Stimuli Response Amplitude');
pP(4) = setParaFields(a{1},'Boolean',0,'relSpeed','Plot Relative Stimuli Response');
pP(5) = setParaFields(a{1},'Boolean',1,'plotErr','Plot Trace SEM Error Bars');
pP(6) = setParaFields(a{2},'Boolean',0,'plotGrid','Plot Metric Bar Graph Gridlines',[],{1,2});
pP(7) = setParaFields(a{1},'Boolean',0,'plotGridT','Plot Stimuli Response Trace Gridlines');
pP(8) = setParaFields(a{2},'Boolean',1,'plotFixedM','Fix Metric Limit To Overall Maximum',[],{1,2});
pP(9) = setParaFields(a{1},'Boolean',0,'plotFixedS','Fix Signal Limit To Overall Maximum');

% if there are multiple stimuli, then set a field list field to select them
if (nDayMx > 1)
    pList2 = cellfun(@(x)(sprintf('Day #%i',x)),num2cell(1:nDayMx)','un',0);
    pP(10) = setParaFields(a{3},'Boolean',0,'grpDayP','Plot All Stimuli Events (Single Type Only)');
    pP(11) = setParaFields(a{3},'List',{1,pList2},'dayInd','Day Index',[],{10,1});    
end

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
[nApp,nSub] = deal(length(snTot.iMov.ok),8);
pF = setFormatFields(nSub);

% initialises the font structs
pF.Title = setFormatFields(setupFontStruct('FontSize',18),'',nSub);
pF.xLabel = setFormatFields(setupFontStruct('FontSize',18),'',nSub);
pF.yLabel = setFormatFields(setupFontStruct('FontSize',18),'',nSub);
pF.Axis = setFormatFields(setupFontStruct('FontSize',12),[],nSub);

% sets x/y labels strings for the 1st plot
pF.Title(1).Font.FontSize = 24;
pF.Axis(1).Font.FontSize = 16;
pF.Title(1).String = 'Stimuli Response';
pF.xLabel(1).String = 'Time Relative To Stimuli Event (min)';
pF.yLabel(1).String = 'Average Speed (mm s^{-1})';

% sets the apparatus names as the titles
pF.Title(2).String = 'GOF Stats';
pF.Title(3).String = 'V_{Amplitude}';
pF.Title(4).String = 'V_{Pre-Stimuli}';
pF.Title(5).String = 'T_{Respond}';
pF.Title(6).String = '\tau_{Activation}';
pF.Title(7).String = '\tau_{Inactivation}';
pF.Title(8).String = '\tau_{Slow}';

% sets the apparatus names as the titles
pF.yLabel(2).String = 'R^2';
pF.yLabel(3).String = 'Speed (mm s^{-1})';
pF.yLabel(4).String = 'Speed (mm s^{-1})';
pF.yLabel(5).String = 'Time (min)';
pF.yLabel(6).String = 'Time (min)';
pF.yLabel(7).String = 'Time (min)';
pF.yLabel(8).String = 'Time (min)';

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Axis(1).String{i} = snTot.iMov.pInfo.gName{i};
end

% sets the font sizes for the inset axes
for i = 2:nSub    
    pF.Title(i).Font.FontSize = 18;
    pF.yLabel(i).Font.FontSize = 12;    
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot,false,false);
[Type1,Type2] = deal(4,2);
[Stats1,Stats2] = deal({'TTestGroup','Tgrp','iDay'},{'GOF','Tgrp','iDay'});
[xDep1,xDep2] = deal({'T','Tgrp','iDay'},{'Tgrp','iDay'});

% sets the independent output variables
oP = addXVarField(oP,'Time','T','Time');
oP = addXVarField(oP,'Time Group','Tgrp','Group');
oP = addXVarField(oP,'Time Group','iDay','Other');

% sets the dependent output variables
oP = addYVarField(oP,'Speed (Abs)','Y',[],Type1,xDep1);
oP = addYVarField(oP,'Speed (Rel)','Y_rel',[],Type1,xDep1);
oP = addYVarField(oP,'Speed (SEM)','Y_sem',[],Type1,xDep1);
oP = addYVarField(oP,'Speed (Fitted)','Y_fit',[],Type1,xDep1);
oP = addYVarField(oP,'Fitted Amp','Yamp',Stats1,[],xDep2);
oP = addYVarField(oP,'Fitted Amp (Mean)','Yamp_mn',[],Type2,xDep2);
oP = addYVarField(oP,'Fitted Amp (SEM)','Yamp_sem',[],Type2,xDep2);
oP = addYVarField(oP,'Raw Amp','YampR',Stats1,[],xDep2);
oP = addYVarField(oP,'Raw Amp (Mean)','YampR_mn',[],Type2,xDep2);
oP = addYVarField(oP,'Raw Amp (SEM)','YampR_sem',[],Type2,xDep2);
oP = addYVarField(oP,'Pre-Stim Spd','YampR',Stats1,[],xDep2);
oP = addYVarField(oP,'Pre-Stim Spd (Mean)','YampR_mn',[],Type2,xDep2);
oP = addYVarField(oP,'Pre-Stim Spd (SEM)','YampR_sem',[],Type2,xDep2);
oP = addYVarField(oP,'Act Tau','kA',Stats1,[],xDep2);
oP = addYVarField(oP,'Act Tau (Mean)','kA_mn',[],Type2,xDep2);
oP = addYVarField(oP,'Act Tau (SEM)','kA_sem',[],Type2,xDep2);
oP = addYVarField(oP,'Inact Tau 1','kI1',Stats1,[],xDep2);
oP = addYVarField(oP,'Inact Tau 1 (Mean)','kI1_mn',[],Type2,xDep2);                        
oP = addYVarField(oP,'Inact Tau 1 (SEM)','kI1_sem',[],Type2,xDep2);                        
oP = addYVarField(oP,'Inact Tau 2','kI2',Stats1,[],xDep2);
oP = addYVarField(oP,'Inact Tau 2 (Mean)','kI2_mn',[],Type2,xDep2);
oP = addYVarField(oP,'Inact Tau 2 (SEM)','kI2_sem',[],Type2,xDep2);
oP = addYVarField(oP,'Max Response Time','Tmax',[],Type2,xDep2);
oP = addYVarField(oP,'Signal GOF','gof',Stats2,[],xDep2);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the parameter struct
cP = retParaStruct(dcObj.pData.cP);

% field retrievals
Tgrp = dcObj.plotD{1}(1).Tgrp;
iAx = dcObj.getSelectAxesIndex;
uStr = {'mm/sec','R^2','mm/sec','mm/sec','min','min','min'};
mStr = {'Speed','GOF','Amplitude','Pre-Stimuli Speed',...
        'Response Time','Activation TC','Inactivation TC'};

% incorporates the fields for the double-exponential field
if cP.useDouble
    mStr{end} = 'Fast Inactivation TC';
    [uStr{end+1},mStr{end+1}] = deal('min','Slow Inactivation TC');
end
    
% sets the common class fields
dcObj.xGrp = Tgrp;
dcObj.yName = mStr{iAx};
dcObj.yUnits = uStr{iAx};
dcObj.combFig = false;
dcObj.grpName = dcObj.pData.appName;

% sets up the axes specfic fields
if iAx == 1
    % case is the speed traces
    dcObj.pType = 'Fitted Trace';  
    dcObj.xName = 'Time';
    dcObj.xUnits = 'sec';
    
else
    % case is the metric bar graphs
    dcObj.pType = 'Bar Graph';
    dcObj.xName = 'Time Bin';
    
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
    % retrieves the parameter struct
    cP = retParaStruct(pData.cP,gPara);
end

% sets the movement calculation type
cP.movType = 'Absolute Speed';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);

% fixed parameters
nGrp = str2double(cP.nGrp);
[devType,chType] = deal([]);
[tBefore,tAfter,iDay,grpDay] = deal(cP.tBefore*60,cP.tAfter*60,[],false);

% setting up of the stimuli signal time vector
T = setStimTimeVector(cP);

% sets the daily time group strings
Tgrp = setTimeGroupStrings(nGrp,cP.Tgrp0);

% retrieves the other calculation parameters (if they exist)
if isfield(cP,'grpDay'); grpDay = cP.grpDay; end
if isfield(cP,'devType'); devType = cP.devType; end
if isfield(cP,'chType'); chType = cP.chType; end

% sets the day index strings (if grouping by day)
if grpDay
    nDayMx = max(detExptDayDuration(snTot));
    iDay = cellfun(@(x)(sprintf('Day #%i',x)),num2cell(1:nDayMx),'un',0);
end

% memory allocation
pData.oP.sepDay = grpDay;
plotD = initPlotValueStruct(snTot,pData,cP,...
                             'T',T','Tgrp',Tgrp,'Hist',[],'iDay',iDay,...
                             'Y',[],'Y_rel',[],'Y_sem',[],'Y_fit',[],...
                             'Yamp_mn',[],'YampR_mn',[],'Y0_mn',[],...
                             'Yamp_sem',[],'YampR_sem',[],'Y0_sem',[],...
                             'kA_mn',[],'kA_sem',[],'kI1_mn',[],'kI1_sem',[],...
                             'kI2_mn',[],'kI2_sem',[],'Tmax',[],'gof',[]);                                             
                                                  
% other memory allocations                         
nDay = 1 + (max(detExptDayDuration(snTot))-1)*grpDay;
indF = cell(nDay,nGrp,nApp);
P = repmat({cell(nDay,nGrp)},nApp,1+(~isempty(snTot(1).Py)));                      
                         
% other initialisations
[hasX,hasY] = deal(false);
[is2D,iOfs] = deal(~isempty(snTot(1).Py),zeros(nApp,1));
aok = cell2mat(cellfun(@(x)(field2cell(...
                    field2cell(x,'iMov',1),'ok')),num2cell(snTot))');
                
% ---------------------------------%
% --- STIMULI TRACE EXTRACTION --- %
% ---------------------------------%

% sets the waitbar offset (is > 0 for more than one
wStr = {'Overall Progress','Fitting Signal Exponentials'};
wOfs = (nExp > 1); wStr = wStr((2-wOfs):end);
if nargin == 5
    % resets the strings
    [h,wOfs] = deal(varargin{1},wOfs+1);
    
    % resets the progressbar fields
    wStr = [h.wStr(1),wStr];
    for i = 2:length(wStr); h.Update(i,wStr{i},0); end
else
    % creates the waitbar figure
    h = ProgBar(wStr,'Stimuli Response Signal Calculation');
end

% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure (if more than one solution file)
    if nExp > 1
        wStrNw = sprintf('%s (Experiment %i of %i)',wStr{wOfs},i,nExp);
        if h.Update(wOfs,wStrNw,i/(1+nExp))
            % if the user cancelled, then exit the function
            [plotD,ok] = deal([],false);
            return
        else
            % resets the other waitbar panels
            h.Update(1+wOfs,'Time Binning Data',0);
        end
    end      
    
    % sets the video/stimuli time stamps into a single vector
    [flyok,Ttot] = deal(snTot(i).iMov.flyok,cell2mat(snTot(i).T));
    [Ts0,~] = getMotorFiringTimes(snTot(i).stimP,devType,chType);
%     Ts = snTot(i).Ts(~cellfun('isempty',snTot(i).Ts));    
    if ~isempty(Ts0)
        % determines the indices of the stimuli events within the total
        % experiment, and determines what time groups that the stimuli
        % events took place in       
        nFlyR = cellfun(@sum,flyok);
        Tsf = Ts0; 
        Ts = num2cell(Tsf);
        indGrp = detTimeGroupIndices(Tsf,...
                        snTot(i).iExpt(1).Timing.T0,nGrp,cP.Tgrp0,grpDay);
                                
        % sets the indices for the pre-stimuli phase
        indV = cellfun(@(x)(find(Ttot>(x-(tBefore+cP.nAvg+1)),1,'first'):...
                            find(Ttot<(x+tAfter),1,'last')),Ts,'un',0);  
        ii = ~cellfun('isempty',indV);               
        
        % sets the total signals for the pre/post stimuli signals
        iTot = cell(size(indV));
        iTot(ii) = cellfun(@(x,y)(roundP((Ttot(x)-y))+...
                    (tBefore+cP.nAvg+1)),indV(ii),Ts(ii),'un',0);                                                                               
                
        % if only considering moving flies, 
        if cP.moveOnly            
            % determines the points for the movement check phase
            indM = cell(size(indV));
            indM(ii) = cellfun(@(x)((find(Ttot<(x-(cP.tMove*60+1)),1,...
                        'last')+1):find(Ttot<x,1,'last')),Ts(ii),'un',0);  
                
            % determines the non-empty regions and allocates memory
            jj = ii & ~cellfun('isempty',indM);
            Pm = repmat({cell(size(indV))},1,2);            
        end                
                
        % calculates the pre/post stimuli velocities for all flies, and
        % bins the values according to their time group bin
        for j = 1:nApp 
            % updates the waitbar figure
            wStrNw = sprintf('%s (Region %i of %i)',...
                             'Time Binning Data',j,nApp);
            if h.Update(1+wOfs,wStrNw,j/(1+nApp))
                % if the user cancelled, then exit the function
                [plotD,ok] = deal([],false);
                return
            end            
            
            if aok(j,i)                 
                % sets the x-location data (if it exists)
                if ~isempty(snTot(i).Px)               
                    % sets the position values for all the signals
                    [Px,hasX] = deal(snTot(i).Px{j}(:,flyok{j}),true);
                    Pxnw = cellfun(@(x,y)(setBinnedSignals(Px(x,:),y,...
                                    length(T)+cP.nAvg)),indV,iTot,'un',0);                      
                    PxG = cellfun(@(x)(Pxnw(x)),indGrp,'un',0);
                    
                    % if only considering moving flies, retrieve the
                    % x-locations of the flies during the moving period
                    if cP.moveOnly
                        Pm{1}(jj) = cellfun(@(x)(Px(x,:)),indM(jj),'un',0);                                                
                    end                                                                 
                end
                
                % sets the y-location data (if it exists)
                if ~isempty(snTot(i).Py)               
                    % sets the position values for all the signals
                    [Py,hasY] = deal(snTot(i).Py{j}(:,flyok{j}),true);
                    Pynw = cellfun(@(x,y)(setBinnedSignals(Py(x,:),y,...
                        length(T)+cP.nAvg)),indV,iTot,'un',0);     
                    PyG = cellfun(@(x)(Pynw(x)),indGrp,'un',0);
                    
                    % if only considering moving flies, retrieve the
                    % x-locations of the flies during the moving period
                    if (cP.moveOnly)
                        Pm{2}(jj) = cellfun(@(x)(Py(x,:)),indM(jj),'un',0);                                                
                    end                                                                                     
                end                
                
                % sets the flags for the flies that are to be kept
                if cP.moveOnly
                    % consider only the moving flies
                    isOK = getPreStimMoving(cP,Pm,jj,nFlyR(j),hasX+2*hasY);
                else
                    % all flies to be considered
                    isOK = repmat({true(1,nFlyR(j))},length(indV),1);                    
                end
                
                % groups the moving fly flags
                isOKG = cellfun(@(x)(isOK(x)),indGrp,'un',0);
                indFT = cellfun(@(x)(cellfun(@(y)(iOfs(j)+find(y)),x,'un',0)),isOKG,'un',0);
                iOfs(j) = iOfs(j) + length(flyok{j});
                
                % converts the sub-cell arrays into numerical arrays
                [iR,iC] = deal(1:size(indFT,1),1:size(indFT,2));
                indF(iR,iC,j) = cellfun(@(x,y)([x,...
                            cell2mat(y(:)')]),indF(iR,iC,j),indFT,'un',0);
                
                % appends the x-locations into the overall array                   
                if (hasX)               
                    PxG = cellfun(@(x,y)(cell2mat(cellfun(@(xx,yy)...
                                (xx(:,yy)),x,y,'un',0)')),PxG,isOKG,'un',0);
                    P{j,1}(iR,iC) = cellfun(@(x,y)([x,y]),P{j,1}(iR,iC),PxG,'un',0); 
                end
                
                % appends the y-locations into the overall array                   
                if (hasY)                    
                    PyG  = cellfun(@(x,y)(cell2mat(cellfun(@(xx,yy)...
                                (xx(:,yy)),x,y,'un',0)')),PyG,isOKG,'un',0);  
                    P{j,2}(iR,iC) = cellfun(@(x,y)([x,y]),P{j,2}(iR,iC),PyG,'un',0);                                            
                end
            end
        end
    end  
end

% loops through each of the apparatus calculating the exponentials
[jj,a,Np] = deal((tBefore+1):length(T),min(1,wOfs),[cP.tBefore,cP.tAfter]*60);
for i = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf('Fitting Exponential Curves (Region %i of %i)',i,nApp);
    if h.Update(1+wOfs,wStrNw,0.5*(1-a)+0.5*(1+a)*(i+1)/(2+nApp))
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    end        
    
    % calculates the velocity trace         
    if any(aok(i,:))
        % calculates the average speed
        [plotD(i).Y,plotD(i).Y_sem] = ...
                        calcSRAvgSpeed(P(i,:),cP,Np,indF(:,:,i),[],is2D);
        plotD(i).Hist = {cellfun(@(x)(size(x,2)),P{i,1})};
                                
        % calculates the initial speed and the relative/absolute stimuli
        % response speed values
        xiT = 1:tBefore;
        plotD(i).Y0_mn = cellfun(@(x)(median(x(xiT,:),1,'omitnan')),...
                                plotD(i).Y,'un',0);        
        plotD(i).Y_rel = cellfun(@(x,y)(x-repmat(y,size(x,1),1)),...
                                plotD(i).Y,plotD(i).Y0_mn,'un',0);

        % fits exponentials to each of the time-binned groups   
        Ydata = cellfun(@(x)(x(jj,:)),plotD(i).Y_rel,'un',0);
        [p,Yfit,gof] = fitSignalExp(plotD(i).T(jj),Ydata,cP);
        plotD(i).gof = cellfun(@(x)(x(:)'),gof,'un',0);
        plotD(i).Y_fit = cellfun(@(x)([zeros(tBefore,nGrp);x]),Yfit,'un',0);
        
        % sets the fitted parameters
        plotD(i) = setSRFittedPara(plotD(i),p,tBefore,nGrp);               
    end
end

% --------------------------------%
% --- HOUSE-KEEPING EXERCISES --- %
% --------------------------------%
    
% closes the waitbar figure
if ~h.Update(1+wOfs,'Stimuli Response Calculations Complete!',1)
    if nargin < 5; h.closeProgBar(); end
end

% ----------------------------------------------------------------------- %
% ---                        PLOTTING FUNCTION                        --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,varargout] = plotFunc(snTot,pData,plotD,ind)

% retrieves the plotting parameter struct
pP = retParaStruct(pData.pP);
sP = retParaStruct(pData.sP);
cP = retParaStruct(pData.cP);
pF = pData.pF;

% sets up any missing fields
if ~isfield(sP,'pT')
    nGrp = str2double(cP.nGrp);
    [sP.pT,sP.pF] = deal(setGroup(1,[nGrp,1]));
end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the indices to be plotted
[iPlotT,iPlotF] = deal(find(sP.pT),find(sP.pF));
if isempty(iPlotT) && isempty(iPlotF); return; end

% sets the number of groups and the parameter struct
p = plotD{1}(sP.pInd);

% if the current data set is empty, then exit the loop
aok = any(cell2mat(cellfun(@(x)(field2cell(...
                field2cell(x,'iMov',1),'ok')),num2cell(snTot))'),2);
if ~aok(sP.pInd); return; end

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% retrieves the day separation popup/checkbox object handles
hPara = findall(0,'tag','figAnalysisPara');
hPopup = findall(hPara,'tag','dayInd');
hCheck = findall(hPara,'tag','grpDayP');
hText = findall(hPara,'String','Day Index: ');

%
if isfield(cP,'grpDay')
    grpDay = cP.grpDay;
else
    grpDay = false;
end

% sets the day index and parameter object properties
if grpDay
    % retrieves the day index flag and group flag    
    if pP.grpDayP
        % if data is grouped by day, but more than one time group is 
        % selected, then output an error message and exit the function
        if (length(iPlotT) > 1) || (length(iPlotF) > 1)
            eStr = ['Error! Unable to plot all daily signal ',...
                    'data for more than one time group.'];
            waitfor(msgbox(eStr,'Analysis Plot Error'))
            return        
        else
            % resets the plotting data struct
            jPlot = unique([iPlotT;iPlotF]);
            p = resetPlotDataStruct(p,true,jPlot);  
            pF.Title(1).String = sprintf('All Days - %s',p.Tgrp0{jPlot});
            
            % resets the plotting indices
            if ~isempty(iPlotT); iPlotT = 1:size(p.Y,2); end
            if ~isempty(iPlotF); iPlotF = 1:size(p.Y,2); end              
        end
    else                
        % otherwise reset the plotting data struct for the selected day
        a = regexp(pP.dayInd,'#','split');
        iDay = str2double(a{end});
        p = resetPlotDataStruct(p,false,iDay);  
        pF.Title(1).String = sprintf('Day %i',iDay);
    end
    
    % updates the object enabled properties
    setObjEnable(hCheck,'on');
    setObjEnable(hPopup,~pP.grpDayP)
    setObjEnable(hText,~pP.grpDayP)
    updateTabEnabledProps(hCheck)
    
else
    % resets the plotting data struct
    p = resetPlotDataStruct(p,false,1);
    
    % if the objects exist, then disable them
    if ~isempty(hPopup)
        set(setObjEnable(hPopup,'off'),'Value',1);
        set(setObjEnable(hCheck,'off'),'Value',false); 
        setObjEnable(hText,'off')
        updateTabEnabledProps(hCheck)
    end
end

% sets the time group value and the x-plot values
nGrp = length(p.Tgrp);
useDouble = ~all(isnan(p.kI2_mn(:)));

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the legend strings
[lStr,colErr] = deal(p.Tgrp,'m');

% sets the main plot ylabel string and offset amount
pF.yLabel(1).String = 'Speed (mm s^{-1})';
if pP.relSpeed
    yOfs = zeros(1,nGrp);    
else
    yOfs = p.Y0_mn;
end

% sets the main plot ylabel string and offset amount
if useDouble
    pF.Title(7).String = '\tau_{Fast}';
end

% sets the title string
grpName = snTot(1).iMov.pInfo.gName{sP.pInd};
pF.Title(1).String = sprintf('%s (%s)',pF.Title(1).String,grpName);

% --------------------------------- %
% --- PLOT AXES INITIALISATIONS --- %
% --------------------------------- %

% sets the fixed dimensions
[nSub,pY,pLim] = deal(6 + useDouble,0.05,0.05*[-1 1]);

% creates the main axes
if pP.incMet
    % creates the sub-plot axes
    hAxIns = cell(nSub,1);
    for i = 1:nSub    
        pPos = calcOuterPos(2,nSub,i+nSub);
        hAxIns{i} = axes('parent',hP,'OuterPosition',pPos,...
                         'Units','Normalized','linewidth',1.5,...
                         'UserData',i+1);
    end
end

% bug error fix - this removes any extraneneous axes from the figure
try
    hAxT = findobj(gcf,'type','axes');
    delete(hAxT(cellfun('isempty',get(hAxT,'UserData'))))
end

% ------------------------------- %
% --- STIMULI RESPONSE TRACES --- %
% ------------------------------- %

% resets the axis position (to shrink if for the axis labels)
if pP.incMet
    oPos = calcOuterPos(2,1,1);
else
    oPos = calcOuterPos(1,1,1);
end

% turns the hold on
hAx = axes('OuterPosition',oPos,'parent',hP,'UserData',1); 
hold(hAx,'on'); 
axis(hAx,'on')
    
% determines the traces/fitted exponentials that are to be plotted
col = getBarColourScheme(nGrp,colErr);
[hPlot,yMin,yMax] = deal(cell(nGrp,1),1e10,-1e10);

% plots the signal traces
for i = 1:length(iPlotT)
    % sets the x/y plot values
    j = iPlotT(i);
    [Tplt,Yplt] = deal(p.T/60,p.Y_rel(:,j)+yOfs(j));
    
    % adds the error bar trace (if selected)
    if pP.plotErr
        plotSignalSEM(Yplt,p.Y_sem(:,j),Tplt,col{j})        
        yMin = min(yMin,min(Yplt-p.Y_sem(:,j)));
        yMax = max(yMax,max(Yplt+p.Y_sem(:,j)));
    else
        [yMin,yMax] = deal(min(yMin,min(Yplt)),max(yMax,max(Yplt)));        
    end    
    
    % creates the plot
    hPlot{j} = plot(Tplt,Yplt,'color',col{j},'UserData',i,'Tag','hRaw');    
end

% plots the exponential traces
for i = 1:length(iPlotF)
    j = iPlotF(i);
    Yplt = p.Y_fit(:,j)+yOfs(j);
    hPlot{j} = plot(p.T/60,Yplt,'color',col{j},...
        'linewidth',2,'UserData',i,'Tag','hFit');
    if pP.relSpeed
        [yMin,yMax] = deal(min(yMin,min(Yplt)),max(yMax,max(Yplt)));        
    else
        [yMin,yMax] = deal(max(0,min(Yplt)),max(yMax,max(Yplt)));        
    end
end

% formats the plot axis
delY = (yMax-yMin);
if abs(delY) < 1e6
    if pP.relSpeed
        set(hAx,'ylim',[roundP(yMin-pY*delY,0.1) roundP(yMax+pY*delY,0.1)]+pLim)
    else
        set(hAx,'ylim',[0 roundP(yMax+pY*delY,0.1)+pLim(2)])
    end
end

% plots the gridlines (if required)
formatPlotAxis(hAx,pF,1);
if pP.plotGridT
    grid(hAx,'on')
end

% sets the axis properties (based on the plot type)
xLim = roundP(p.T([1 end])/60);
set(hAx,'xlim',xLim,'linewidth',1.5,'box','on'); 
if pP.relSpeed; plot(hAx,xLim,[0 0],'r--','linewidth',1.5); end  

% retrieves the panel object handle
yLim = getCurrentAxesProp('ylim');
plot(hAx,[0 0],yLim,'r--','linewidth',1.5)

% if not plotting the metrics, then exit the function
if ~pP.incMet                      
    % creates the legend object using the formatting data struct    
    ii = ~cellfun('isempty',hPlot);
    if any(ii) && (nGrp > 1)
        % creates the legend
        pF.Legend.String = cellfun(@(x,y)(sprintf('%s (N = %i)',x,y)),...
                        lStr(ii),num2cell(p.Hist(ii)),'un',0);

        % creates the new legend    
        hLg = createLegendObj(hPlot(ii),pF.Legend,1,0);
    end        
    
    % exits the function
    return
end

% ---------------------------------------- %
% --- EXPONENTIAL FIT PARAMETER INSETS --- %
% ---------------------------------------- %

% initialisations
[xLim,Wax,Ymx,maxAxR] = deal(get(hAx,'xlim'),0.97,zeros(1,nSub),100);

% sets the upper/lower limits
pD = plotD{1}(~cellfun('isempty',field2cell(plotD{1},'Y')));
if (pP.plotFixedM)
    % amplitude maximum limits
    if (pP.plotRaw)
        [YampMn,YampS] = field2cell(pD,{'YampR_mn','YampR_sem'});            
    else
        [YampMn,YampS] = field2cell(pD,{'Yamp_mn','Yamp_sem'});
    end

    % sets the amplitude maximum limits
    Ymx(2) = max(cellfun(@(xx,yy)(max(cellfun...
                    (@(x,y)(detOverallLimit(x+y)),xx,yy))),YampMn,YampS));    
    
    % sets the pre-stimuli average speed maximum limits
    [Y0Mn,Y0S] = field2cell(pD,{'Y0_mn','Y0_sem'});
    Ymx(3) = max(cellfun(@(xx,yy)(max(cellfun...
                    (@(x,y)(detOverallLimit(x+y)),xx,yy))),Y0Mn,Y0S));                       
end

% sets the new left location of the main trace plot
[pLbl,pAx] = deal(get(get(hAx,'yLabel'),'Extent'),get(hAx,'position'));
pAx(1) = pAx(3)*(xLim(1)-pLbl(1))/diff(xLim);

% plots the inset metrics (if selected)
if pP.incMet
    % sets the parameter strings
    pStr = {'gof','Yamp_mn','Y0_mn','Tmax','kA','kI1','kI2'};
    if pP.plotRaw; pStr{2} = 'YampR_mn'; end
    pStr = pStr(1:nSub);    
    
    % sets the axis properties for each of the inset plots
    for i = 1:length(pStr)   
        % set the hold on for the subplot
        hold(hAxIns{i},'on');         
        
        % retrieves the stimuli response values
        [Ynw,YSEM,YmxNw] = getSRValues(p,pD,pP,pStr{i});        
        if ~isnan(YmxNw); Ymx(i) = YmxNw; end
        
        % creates the bars for each of the plot values
        for j = 1:nGrp
            bar(hAxIns{i},j,Ynw(j),'facecolor',col{j},'tag','hBar','UserData',j)
        end
        
        % sets the axis labels
        set(hAxIns{i},'xticklabel',[],'xtick',1:nGrp,'linewidth',1.5,...
                      'xlim',[1 nGrp] + 0.5*[-1 1],'Box','on')        

        % sets the axis limits
        if strcmp(pStr{i},'gof')
            set(hAxIns{i},'ylim',[0 1])        
        else            
            if pP.plotFixedM
                resetYAxisScale(hAxIns{i},Ynw(:),YmxNw);  
            elseif isempty(YSEM)
                resetYAxisScale(hAxIns{i},Ynw(:));
            else
                resetYAxisScale(hAxIns{i},Ynw(:)+YSEM(:));
            end
        end
        
        % formats the plot axis
        set(gcf,'CurrentAxes',hAxIns{i});
        pF.xLabel(i+1).String = pF.Axis(i+1).String;
        formatPlotAxis(hAxIns{i},pF,1+i);      
        set(hAxIns{i},'xtick',nGrp/2+0.5)
        plot(hAxIns{i},get(hAxIns{i},'xlim'),[0 0],'k','linewidth',2)
        
        % if there is a SEM signal to add, then add it...
        if ~isempty(YSEM)
            addBarError(hAxIns{i},1:nGrp,Ynw,YSEM,colErr,2);
        end        
        
        % turns the grid on (if specified)
        if pP.plotGrid
            set(hAxIns{i},'ygrid','on')
        end
    end
    
    % resets the GOF stat sub-plot so it matches the others
    axP = get(hAxIns{2},'position');
    resetObjPos(hAxIns{1},'Height',axP(4));
    
    % creates the legend object
    if nGrp > 1
        % resets the legend strings
        pF.Legend.String = cellfun(@(x,y)(sprintf('%s (N = %i)',x,y)),...
                        lStr(:),num2cell(p.Hist(:)),'un',0);    
           
        % resets the bar objects
        hBar = cell(nGrp,1);
        hBar0 = findall(get(hAxIns{1},'Parent'),'tag','hBar');

        for j = 1:nGrp
            hBarNw = findall(hBar0,'UserData',j);
            hBar{j} = hBarNw(1);
        end      

        % creates the legend object
        hLg = createLegendObj(hBar,pF.Legend,1,0);          
        
        % creates the legend object       
        lgP = resetVertLegendWidth(hLg);
        lgP(1:2) = [(1-lgP(3)),(0.75-lgP(4)/2)];
        set(hLg,'position',lgP);      
        
        % resets the width of the trace plot axis
        lgPNw = get(hLg,'position');
        pAx(3) = lgPNw(1) - pAx(1);
    else
        % sets the width to the full extent        
        pAx(3) = Wax - pAx(1);
    end    
else
    % sets the width to the full extent
    pAx(3) = Wax - pAx(1);    
end

% resets the trace plot position
set(hAx,'position',pAx)

% formats the plot axis to the fitted response
if any(sP.pT) || any(sP.pF) 
    % sets the axis properties (based on the plot type) 
    if pP.plotFixedS
        [Y,Ysem,Y0_mn] = field2cell(plotD{1},{'Y','Y_sem','Y0_mn'});
        if ~pP.relSpeed
            Y = cellfun(@(xx,yy)(cellfun(@(x,y)(x+repmat...
                        (y,size(x,1),1)),xx,yy,'un',0)),Y,Y0_mn,'un',0);
        end
        
        % determines the overall min/max y-axis limits
        yLim0 = get(hAx,'ylim');
        Ymx = max(cellfun(@(xx,yy)(cellfun(@(x,y)(max(x(:)+y(:))),xx,yy)),Y,Ysem));        
        Ymn = min(cellfun(@(xx,yy)(cellfun(@(x,y)(min(x(:)+y(:))),xx,yy)),Y,Ysem));        
        yLim = [min(Ymn,yLim0(1)) max(Ymx,yLim0(2))];
    else        
        % otherwise, use the current axis limits
        yLim = get(hAx,'ylim');
    end

    % plots the stimuli marker line
    if pP.plotErr     
        plot(hAx,[0 0],yLim,'r--','linewidth',1.5,'HitTest','off')        
        set(hAx,'linewidth',1.5,'ylim',yLim)
    else
        plot(hAx,[0 0],get(hAx,'ylim'),'r--','linewidth',2,'HitTest','off')
    end
    
    % plots a line for zero speed (relative speed only)
    if pP.relSpeed
        xLim = getCurrentAxesProp('xlim');
        plot(hAx,xLim,[0 0],'r--','linewidth',1.5,'HitTest','off')
    end
end

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% initialisations
[cP,sepDay] = deal(retParaStruct(pData.cP),false);

% retrieves the day separation flag (if available)
if (isfield(cP,'grpDay'))
    sepDay = cP.grpDay;
end

% removes the day index dependency if not separating by day
if (~sepDay)
    % determines which output variables have a day index dependency
    xDep = field2cell(pData.oP.yVar,'xDep');
    ii = find(cellfun(@(x)(any(strcmp(x,'iDay'))),xDep));
    
    % removes the day index dependency
    for i = 1:length(ii)
        pData.oP.yVar(ii(i)).xDep = xDep{ii(i)}(~strcmp(xDep{ii(i)},'iDay'));
    end
    
    % determines which output variables have a day index dependency
    Stats = field2cell(pData.oP.yVar,'Stats');
    jj = find(cellfun(@(x)(any(strcmp(x,'iDay'))),Stats));
    
    % removes the day index dependency
    for i = 1:length(jj)
        pData.oP.yVar(jj(i)).Stats = Stats{jj(i)}(~strcmp(Stats{jj(i)},'iDay'));
        if (strcmp(pData.oP.yVar(jj(i)).Stats{1},'TTestGroup'))
            pData.oP.yVar(jj(i)).Stats{1} = 'TTest';    
        end
    end    
end

% resets the group separation flag
pData.oP.sepGrp = sepDay;

% ----------------------------------------------------------------------- %
% ---                         OTHER FUNCTIONS                         --- %
% ----------------------------------------------------------------------- %

% --- determines if the flies have moved during the pre-stimuli phase
function isMove = getPreStimMoving(cP,Pm,jj,nFlyR,Type)

% memory allocation
isMove = repmat({true(1,nFlyR)},length(Pm{1}),1);

% determines the fly movement range (based on the type)
switch (Type)
    case (1) % case is 1D (x-locations only)
        D = cellfun(@(x)(range(x,1)),Pm{1}(jj),'un',0);
    case (2) % case is 1D (y-locations only)
        D = cellfun(@(x)(range(x,1)),Pm{2}(jj),'un',0);
    case (3) % case is 2D
        D = cellfun(@(x,y)(max(range(x,1),range(y,1))),Pm{1}(jj),Pm{2}(jj),'un',0);
end

% determines if the flies have moved through the required distance
isMove(jj) = cellfun(@(x)(x > cP.dMove),D,'un',0);

% --- 
function p = resetPlotDataStruct(p,grpDay,ind)

% sets the fields not to be altered
pStr = {'T','Tgrp','iDay'};

% retrieves the plotting data struct field names (removes the fields above)
fName = fieldnames(p);
fName = fName(cellfun(@(x)(~any(strcmp(x,pStr))),fName));

% resets the plotting data struct values
if (grpDay)
    % case is plotting output data for all days
    for j = 1:length(fName)
        A = eval(sprintf('p.%s',fName{j}));
        Ynw = cell2mat(cellfun(@(x)(x(:,ind)),A,'un',0)');
        eval(sprintf('p.%s = Ynw;',fName{j}));
    end
    
    % resets the time group strings
    p.Tgrp0 = p.Tgrp;
    p.Tgrp = cellfun(@(x)(sprintf('Day #%i',x)),num2cell(1:size(p.Y,2))','un',0);
else
    % case is plotting output data for a single day
    for j = 1:length(fName)
        if (strcmp(fName{j},'Hist'))
            p.Hist = p.Hist{1}(ind,:);
        else
            eval(sprintf('p.%s = p.%s{ind};',fName{j},fName{j}));
        end
    end
end