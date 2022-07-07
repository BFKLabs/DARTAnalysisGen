% --- calculation of the fly activity exponential response with respect to 
%     delivered stimuli (short experiment only)
function pData = StimResponseShort(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Stimuli Response Curve Fitting (Short)';
pData.Type = {'Pop','Multi'};
pData.fType = [2 2 2 1];
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
    [pData.hasSR,pData.hasRS] = deal(true,false);
    pData.sP = initSpecialPara(snTotL,pData,'appNameS',1); 
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
rI.Dur = 'Short';
rI.Shape = 'None';
rI.Stim = 'Motor';
rI.Spec = 'None';

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% determines if the stimuli protocols for each experiment are identical
canGrpStim = hasEqualStimProtocol(snTot);

% initialises the parameter struct
nPara = 9 + canGrpStim;                       
cP = setParaFields(nPara);

% tab list headings
a = {'1 - General','2 - Exponential Fit'};

% sets the parameter list for the drop-down list
pList = {'Mean','Median'};

% sets the parameter fields
cP(1) = setParaFields(a{1},'List',{1,pList},'cType','Speed Calculation Type');
cP(2) = setParaFields(a{1},'Number',2,'tBefore','Time Before Stimuli (min)',[1 15 true]);
cP(3) = setParaFields(a{1},'Number',15,'tAfter','Time After Stimuli (min)',[1 45 true]);
cP(4) = setParaFields(a{1},'Boolean',0,'moveOnly','Ignore Non-Reactive Flies For Calculations');
cP(5) = setParaFields(a{1},'Number',5,'tMove','Non-Reactivity Duration (min)',[1 60 true],{4,2});
cP(6) = setParaFields(a{2},'Boolean',0,'useDouble','Use Double Exponential Inactivity Fit');
cP(7) = setParaFields(a{2},'Boolean',0,'useToff','Use Exponential Time Offset');
cP(8) = setParaFields(a{2},'Boolean',0,'ignoreWeak','Ignore Weak Signal Response');
cP(9) = setParaFields(a{2},'Number',0.001,'pWeak','Peak/Signal Length Ratio', [0 0.1 true],{8,2});

% sets the tool-tip strings
cP(1).TTstr = 'The speed calculation type (mean or median)';
cP(2).TTstr = 'The signal duration before the stimuli event';
cP(3).TTstr = 'The signal duration after the stimuli event';
cP(4).TTstr = 'Duration for which a fly is considered to be non-reactive';
cP(5).TTstr = 'Removes non-reactive flies from the stimuli response signal';
cP(6).TTstr = 'Uses double exponential fit instead of single exponential';
cP(7).TTstr = 'Includes a time offset for the exponential fit';
cP(8).TTstr = 'Ignore signals that are too weak to be analysed';
cP(9).TTstr = 'Peak/Signal length ratio beneath which the signal is considered too weak';

% sets the parameter fields/tt-strings if data can be grouped by stimuli 
if (canGrpStim)
    % sets the parameter fields
    cP(10) = setParaFields(a{1},'Boolean',0,'grpStim','Separate Data By Stimuli Events',[],{{{'SRS'}},{1}});
        
    % sets the tool-tip strings
    cP(10).TTstr = 'Flag indicating whether the response is to be grouped by stimuli events';    
end

% adds the unique motor parameters
cP = addUniqueMotorPara(cP,snTot);

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% % determines if the stimuli protocols for each experiment are identical
% [canGrpStim,nStim] = hasEqualStimProtocol(snTot);

% initialises the parameter struct
nPara = 7;                      
pP = setParaFields(nPara);

% tab list headings
a = {'1 - Stimuli Trace','2 - Fit Metrics','3 - Trace Grouping'};

% sets the parameter list for the drop-down list
pList = {'R-Squared','Adjusted R-Squared'};

% sets the parameter fields
pP(1) = setParaFields(a{2},'Boolean',1,'incMet','Include Stimuli Response Metric Subplots');
pP(2) = setParaFields(a{2},'List',{1,pList},'gofType','Goodness-of-fit Type',[],{1,2});
pP(3) = setParaFields(a{1},'Boolean',1,'plotRaw','Plot Raw Trace Stimuli Response Amplitude');
pP(4) = setParaFields(a{1},'Boolean',0,'relSpeed','Plot Relative Stimuli Response');
pP(5) = setParaFields(a{1},'Boolean',1,'plotErr','Plot Trace SEM Error Bars');
pP(6) = setParaFields(a{2},'Boolean',0,'plotGrid','Plot Metric Bar Graph Gridlines',[],{1,2});
pP(7) = setParaFields(a{1},'Boolean',0,'plotGridT','Plot Stimuli Response Trace Gridlines');

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
pF.Legend = setFormatFields(setupFontStruct('FontSize',12),[],1);

% sets x/y labels strings for the 1st plot
pF.Title(1).Font.FontSize = 24;
pF.Axis(1).Font.FontSize = 16;
pF.Title(1).String = 'Stimuli Response Trace';
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
function oP = initOutputPara(snTot,cP)

% initialisations
oP = setupOutputParaStruct(snTot,false,false);
[Stats1,Stats2] = deal({'CompSumm','iStimS'},{'GOF','iStimS'});
[Type1,Type2,xDep1,xDep2] = deal(4,2,{'T','iStimS'},{'iStimS'});

% sets the independent output variables
oP = addXVarField(oP,'Time','T','Time');
oP = addXVarField(oP,'Time Group','iStimS','Other');

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
oP = addYVarField(oP,'Pre-Stim Spd','Y0',Stats1,[],xDep2);
oP = addYVarField(oP,'Pre-Stim Spd (Mean)','Y0_mn',[],Type2,xDep2);
oP = addYVarField(oP,'Pre-Stim Spd (SEM)','Y0_sem',[],Type2,xDep2);
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

% initialises the output parameter struct
pData.oP = initOutputPara(snTot,cP);

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% memory allocation
is2D = ~isempty(snTot(1).Py);
[tBefore,tAfter] = deal(cP.tBefore*60,cP.tAfter*60);
T = (-tBefore:tAfter)';
nExp = length(snTot);

% determines the stimuli event count (based on whether they are grouped)
if (isfield(cP,'grpStim'))
    grpStim = cP.grpStim;
else
    grpStim = false;
end
    
% sets the stimuli count/indices based on the parameter selection
if (grpStim)
    % stimuli separation has been selected
    [~,nStim] = hasEqualStimProtocol(snTot);
    iStim = (1:nStim);
    iStimS = cellfun(@(x)(sprintf('Stim #%i',x)),num2cell(iStim),'un',0)';
else
    % stimuli separation not selected/possible
    [nStim,iStim,iStimS] = deal(1,[],[]);
end

% retrieves the other calculation parameters (if they exist)
[dT,chT] = deal([]);
if isfield(cP,'devType'); dT = cP.devType; end
if isfield(cP,'chType'); chT = cP.chType; end

% memory allocation
aok = cell2mat(cellfun(@(x)(field2cell(...
                    field2cell(x,'iMov',1),'ok')),num2cell(snTot))');
nApp = length(snTot(1).iMov.ok);
plotD = initPlotValueStruct(snTot,pData,cP,'gType',grpStim,...
                         'T',T,'iStim',iStim,'iStimS',iStimS,'Hist',[],...
                         'Y',[],'Y_rel',[],'Y_sem',[],'Y_fit',[],...
                         'Yamp_mn',[],'YampR_mn',[],'Y0_mn',[],...
                         'Yamp_sem',[],'YampR_sem',[],'Y0_sem',[],...
                         'kA_mn',[],'kA_sem',[],'kI1_mn',[],'kI1_sem',[],...
                         'kI2_mn',[],'kI2_sem',[],'Tmax',[],'gof',[]);                                             
                       
% parameters
[hasX,hasY] = deal(false);
[iOfs,indF,indS] = deal(zeros(nApp,1),cell(1,nStim,nApp),cell(nApp,1));
P = repmat({cell(1,nStim)},nApp,1+(~isempty(snTot(1).Py)));

% other initialisations
cP.nGrp = num2str(nStim);

% retrieves the stimuli start/finish times
[stimP,sTrain] = field2cell(snTot,{'stimP','sTrainEx'});
Ts = cellfun(@(x)(getMotorFiringTimes(x,dT,chT)),stimP,'un',0);
% [Ts,Tf] = getStimTimes(stimP,sTrain,'Motor');
% [Ts,Tf] = deal(Ts{1},Tf{1});

% ---------------------------------%
% --- STIMULI TRACE EXTRACTION --- %
% ---------------------------------%

% sets the waitbar offset (is > 0 for more than one
wStr = {'Overall Progress','Stimuli Signal Binning'};
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
            [plotD,ok] = deal([],false);
            return
        end
    end
        
    % sets the video/stimuli time stamps into a single vector
    [flyok,Ttot] = deal(snTot(i).iMov.flyok,cell2mat(snTot(i).T));   
%     Ts = snTot(i).Ts(~cellfun(@isempty,snTot(i).Ts));    
    if ~isempty(Ts{i})
        % sets the indices for the pre/post stimuli phases
        nFlyR = cellfun(@sum,flyok);
        indV = arrayfun(@(x)(find(Ttot>(x-(tBefore+cP.nAvg+1)),1,...
                    'first'):find(Ttot<(x+tAfter),1,'last')),Ts{i},'un',0);               

        % sets the total signals for the pre/post stimuli signals
        iTot = cellfun(@(x,y)(roundP((Ttot(x)-y))+...
                        (tBefore+cP.nAvg+1)),indV,num2cell(Ts{i}),'un',0);
                        
        % if only considering moving flies, 
        if cP.moveOnly               
            % determines the points for the movement check phase
            indM = arrayfun(@(x)((find(Ttot<(x-(cP.tMove*60+1)),1,...
                        'last')+1):find(Ttot<x,1,'last')),Ts{i},'un',0);  
                
            % determines the non-empty regions and allocates memory
            jj = ~cellfun(@isempty,indM);
            Pm = repmat({cell(size(indV))},1,2);            
        end                      
                    
        % calculates the pre/post stimuli velocities for all flies, and
        % bins the values according to their time group bin
        for j = 1:nApp  
            % updates the waitbar figure
            wStrNw = sprintf('Time Binning Data (Region %i of %i)',j,nApp);
            if h.Update(1+wOfs,wStrNw,j/nApp)
                % if the user cancelled, then exit the function
                [plotD,ok] = deal([],false);
                return
            end               
            
            if aok(j,i) 
                % sets the x-location data (if it exists)                
                if ~isempty(snTot(i).Px{j})               
                    % sets the position values for all the signals
                    [Px,hasX] = deal(snTot(i).Px{j}(:,flyok{j}),true);
                    Pxnw = cellfun(@(x,y)(setBinnedSignals(Px(x,:),y,...
                        length(T)+cP.nAvg)),indV,iTot,'un',0);

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
                
                %
                if nStim > 1
                    indST = cellfun(@(x,y)(y*ones(1,sum(x))),...
                                    isOK,num2cell(1:length(Ts))','un',0);
                    indS{j} = [indS{j},cell2mat(indST(:)')];
                end
                                
                % sets the fly indices and 
                indFT = cellfun(@(x)(find(x)+iOfs(j)),isOK(:)','un',0);
                iOfs(j) = iOfs(j) + nFlyR(j);
                
                % converts the sub-cell arrays into numerical arrays
                if grpStim
                    indF(:,:,j) = cellfun(@(x,y)([x,y(:)']),...
                                            indF(:,:,j),indFT,'un',0);                
                else
                    indF{1,1,j} = [indF{1,1,j},cell2mat(indFT)];
                end
                                        
                % appends the x-locations into the overall array                   
                if hasX
                    PxG = cellfun(@(x,y)(x(:,y)),Pxnw,isOK,'un',0)';
                    if ~grpStim; PxG = {cell2mat(PxG)}; end
                    P{j,1} = cellfun(@(x,y)([x,y]),P{j,1},PxG,'un',0); 
                end
                
                % appends the y-locations into the overall array    
                if hasY                
                    PyG = cellfun(@(x,y)(x(:,y)),Pynw,isOK,'un',0)';
                    if ~grpStim; PyG = {cell2mat(PyG)}; end
                    P{j,2} = cellfun(@(x,y)([x,y]),P{j,2},PyG,'un',0); 
                end
            end
        end
    end  
end

% loops through each of the apparatus calculating the exponentials
Np = [cP.tBefore,cP.tAfter]*60;
[jj,a] = deal((tBefore+1):length(T),max(1,wOfs));
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
        % calculates the average speed from the signals
        [plotD(i).Y,plotD(i).Y_sem] = ...
                        calcSRAvgSpeed(P(i,:),cP,Np,indF(:,:,i),[],is2D);
        plotD(i).Hist = {cellfun(@(x)(size(x,2)),P{i,1})};
                    
        if ~iscell(plotD(i).Y)
            plotD(i).Y = {plotD(i).Y};
            plotD(i).Y_sem = {plotD(i).Y_sem};
        end
        
        % calculates the initial speed and the relative/absolute stimuli
        % response speed values
        xiT = 1:tBefore;
        plotD(i).Y0_mn = cellfun(@(x)(median(x(xiT,:),1,'omitnan')),...
                            plotD(i).Y,'un',0);        
        plotD(i).Y_rel = cellfun(@(x,y)(x-repmat(y,size(x,1),1)),...
                            plotD(i).Y,plotD(i).Y0_mn,'un',0);  
        
        % fits exponentials to each of the time-binned groups   
        Ydata = plotD(i).Y_rel{1}(jj,:);
        [p,Yfit,gof] = fitSignalExp(plotD(i).T(jj),Ydata,cP);
        plotD(i).gof = {gof'};
        plotD(i).Y_fit = {[zeros(tBefore,nStim);Yfit]};
        
        % sets the fitted parameters
        plotD(i) = setSRFittedPara(plotD(i),{p},tBefore,nStim);  
    else
        % creates an empty GOF Stats array
        plotD(i).gof = repmat(setEmptyGOF,[1 nStim]);
    end
end

% --------------------------------%
% --- HOUSE-KEEPING EXERCISES --- %
% --------------------------------%

% closes the waitbar figure
if ~h.Update(wOfs+1,'Stimuli Response Calculations Complete!',1)
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

% determines the stimuli event count (based on whether they are grouped)
if (isfield(cP,'grpStim'))
    grpStim = cP.grpStim;
else
    grpStim = false;
end

% if the stimuli grouping checkbox flag does not match the calculated data,
% then output an error to screen and exit the function
if (grpStim ~= plotD{1}(1).gType)
    eStr = {'The stimuli separation flag does not match the calculated data.';'';...
            'Either toggle the separation flag checkbox or recalculate the data.'};
    waitfor(msgbox(eStr,'Invalid Parameter Selection','modal'))
    return
end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the indices to be plotted
[iPlotT,iPlotF] = deal(find(sP.pT),find(sP.pF));
if (isempty(iPlotT) && isempty(iPlotF)); return; end

% sets the plot data struct
p = resetPlotDataStruct(plotD{1},grpStim,sP.pInd);
nApp = size(p.Y,2);

% sets the group indices (based on the stimuli grouping type)
if (grpStim)
    % stimuli events separation selected
    iApp = sP.pInd;
else
    % stimuli events separation not selected/possible
    iApp = unique([iPlotT(:);iPlotF(:)]);
end

% sets the apparatus count
[colErr,useDouble] = deal('m',~all(isnan(p.kI2_mn(:))));

% if the current data set is empty, then exit the function
aok = any(cell2mat(cellfun(@(x)(field2cell(...
                field2cell(x,'iMov',1),'ok')),num2cell(snTot))'),2);
if (~any(aok(unique(iApp)))); return; end    

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the legend strings
if grpStim
    grpStr = snTot(1).iMov.pInfo.gName{iApp};
    lStr = arrayfun(@(x)(sprintf('Stimuli #%i',x)),(1:nApp)','un',0);    
    pF.Title(1).String = sprintf('Stimuli Response (%s)',grpStr);
else
    lStr = snTot(1).iMov.pInfo.gName';    
end

% sets the main plot ylabel string and offset amount
pF.yLabel(1).String = 'Speed (mm s^{-1})';
if (pP.relSpeed)
    yOfs = zeros(1,nApp);    
else
    yOfs = p.Y0_mn;
end    
    
% sets the main plot ylabel string and offset amount
if (useDouble)
    pF.Title(7).String = '\tau_{Fast}';
end

% if the activation time constants are small, then convert to seconds
if (all(p.kA_mn(~isnan(p.kA_mn(iApp))) < 5e-3))
    p.kA_mn = p.kA_mn*60;
    pF.yLabel(6).String = 'Time (sec)';        
end
    
% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% --------------------------------- %
% --- PLOT AXES INITIALISATIONS --- %
% --------------------------------- %

% sets the fixed dimensions
nSub = 6 + useDouble;

% creates the main axes
if (pP.incMet)
    % creates the sub-plot axes
    hAxIns = cell(nSub,1);
    for i = 1:nSub    
        outerPos = calcOuterPos(2,nSub,i+nSub);
        hAxIns{i} = axes('parent',hP,'Units','Normalized','linewidth',...
                        1.5,'UserData',i,'OuterPosition',outerPos);
    end
end

% bug error fix - this removes any extraneneous axes from the figure
try
    hAxT = findobj(gcf,'type','axes');
    delete(hAxT(cellfun(@isempty,get(hAxT,'UserData'))))
end

% ------------------------------- %
% --- STIMULI RESPONSE TRACES --- %
% ------------------------------- %

% resets the axis position (to shrink if for the axis labels)
if (pP.incMet)
    oPos = calcOuterPos(2,1,1);
else
    oPos = calcOuterPos(1,1,1);
end

% turns the hold on
hAx = axes('OuterPosition',oPos,'parent',hP); 
hold(hAx,'on'); axis(hAx,'on')

% determines the traces/fitted exponentials that are to be plotted
col = getBarColourScheme(nApp,colErr);
hPlot = cell(nApp,1);
    
% plots the raw signal traces
for i = 1:length(iPlotT)
    % sets the x/y plot values
    j = iPlotT(i);
    [Tplt,Yplt] = deal(p.T/60,p.Y_rel(:,j)+yOfs(j));
    
    % adds the error bar trace (if selected)
    if (pP.plotErr)
        plotSignalSEM(Yplt,p.Y_sem(:,j),Tplt,col{j})        
    end    
    
    % creates the plot
    hPlot{j} = plot(Tplt,Yplt,'color',col{j});    
end

% plots the fitted exponential traces
for i = 1:length(iPlotF)
    % sets the x/y plot values
    j = iPlotF(i);
    Yplt = p.Y_fit(:,j)+yOfs(j);
    
    % creates the plot
    hPlot{j} = plot(p.T/60,Yplt,'color',col{j},'linewidth',2);
end

% plots the gridlines (if required)
formatPlotAxis(hAx,pF,1);
if (pP.plotGridT)
    grid(hAx,'on')
end

% sets the axis properties (based on the plot type)
xLim = roundP(p.T([1 end])/60);
set(hAx,'xlim',xLim,'linewidth',1.5,'box','on'); 
if (pP.relSpeed); plot(hAx,xLim,[0 0],'r--','linewidth',1.5); end 

% plots the stimuli marker (ensures minimum is at 0)
yL = get(hAx,'ylim'); yL(1) = min(yL(1),0); 
plot(hAx,[0 0],yL,'r--','linewidth',1.5,'tag','hStim')

% if not plotting the metrics, then exit the function
if (~pP.incMet)                           
    % creates the legend object using the formatting data struct    
    ii = ~cellfun(@isempty,hPlot);
    if (any(ii) && (nApp > 1))
        % sets the legend strings
        pF.Legend.String = cellfun(@(x,y)(sprintf('%s (N = %i)',x,y)),...
                                lStr(ii),num2cell(p.Hist(ii))','un',0);

        % creates the new legend    
        hLg = createLegendObj(hPlot(ii),pF.Legend,1,0);

        % creates the legend object 
        [lgP,dX] = deal(resetVertLegendWidth(hLg),0.025);
        lgP(1:2) = [(1-(lgP(3)+dX)),(0.5-lgP(4)/2)];
        set(hLg,'position',lgP);      

        % resets the width of the trace plot axis
        pAx = get(hAx,'position');
        lgPNw = get(hLg,'position');
        resetObjPos(hAx,'width',lgPNw(1) - (pAx(1) + dX/2));
    end        
    
    % exits the function
    return
end

% ---------------------------------------- %
% --- EXPONENTIAL FIT PARAMETER INSETS --- %
% ---------------------------------------- %
 
% initialisations
[xLim,Wax,maxAxR] = deal(get(hAx,'xlim'),0.97,100);

% sets the new left location of the main trace plot
pD = plotD{1}(~cellfun(@isempty,field2cell(plotD{1},'Y')));
[pLbl,pAx] = deal(get(get(hAx,'yLabel'),'Extent'),get(hAx,'position'));
pAx(1) = pAx(3)*(xLim(1)-pLbl(1))/diff(xLim);

% plots the inset metrics (if selected)
if (pP.incMet)
    % sets the parameter strings
    pStr = {'gof','Yamp_mn','Y0_mn','Tmax','kA','kI1','kI2'};
    if (pP.plotRaw); pStr{2} = 'YampR_mn'; end
    pStr = pStr(1:nSub);    
    
    % sets the axis properties for each of the inset plots
    for i = 1:length(pStr)
        % set the hold on for the subplot
        hold(hAxIns{i},'on');         
        
        % sets the SEM signal 
        [Ynw,YSEM,~] = getSRValues(p,pD,pP,pStr{i});
        
        % creates the bars for each of the plot values
        for j = 1:nApp
            bar(hAxIns{i},j,Ynw(j),'facecolor',col{j},'tag','hBar','UserData',j)
        end            

        % sets the axis labels
        set(hAxIns{i},'xticklabel',[],'xtick',1:nApp,'linewidth',1.5,...
                      'xlim',[1 nApp] + 0.5*[-1.05 1.05])        
                
        % formats the plot axis
        pF.xLabel(i+1).String = pF.Axis(i+1).String;
        formatPlotAxis(hAxIns{i},pF,1+i);      
        set(hAxIns{i},'xtick',nApp/2+0.5)
        plot(hAxIns{i},get(hAxIns{i},'xlim'),[0 0],'k','linewidth',2)
        
        % if there is a SEM signal to add, then add it...
        if (strcmp(pStr{i},'gof'))
            set(hAxIns{i},'ylim',[0 1]) 
        elseif (~isempty(YSEM))
            addBarError(hAxIns{i},1:nApp,Ynw,YSEM,colErr,2);
            if (all(isnan(YSEM))); YSEM(:) = 0; end
            resetYAxisScale(hAxIns{i},Ynw(:)+YSEM(:));            
        else
            resetYAxisScale(hAxIns{i},Ynw(:));
        end        
        
        % turns the grid on (if specified)
        if (pP.plotGrid)
            set(hAxIns{i},'ygrid','on')
        end
    end    
    
    % creates the legend object
    if (nApp > 1)
        % sets the legend strings
        pF.Legend.String = cellfun(@(x,y)(sprintf('%s (N = %i)',x,y)),...
                            lStr(:),num2cell(p.Hist(:)),'un',0);                     
        hBar = cell(nApp,1);
        hBar0 = findall(get(hAxIns{1},'Parent'),'tag','hBar');

        for j = 1:nApp
            hBarNw = findall(hBar0,'UserData',j);
            hBar{j} = hBarNw(1);
        end
        
        % creates the legend object
        hLg = createLegendObj(hBar,pF.Legend,1,0);          
        
        % creates the legend object 
        [lgP,dX] = deal(resetVertLegendWidth(hLg),0.025);
        lgP(1:2) = [(1-(lgP(3)+dX)),(0.75-lgP(4)/2)];
        set(hLg,'position',lgP);      
        
        % resets the width of the trace plot axis
        lgPNw = get(hLg,'position');
        pAx(3) = lgPNw(1) - (pAx(1) + dX/2);
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

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% initialisations
[cP,sepStim] = deal(retParaStruct(pData.cP),false);

% retrieves the day separation flag (if available)
if (isfield(cP,'grpStim'))
    sepStim = cP.grpStim;
end

% removes the day index dependency if not separating by day
if (~sepStim)
    % determines which output variables have a day index dependency
    xDep = field2cell(pData.oP.yVar,'xDep');
    ii = find(cellfun(@(x)(any(strcmp(x,'iStimS'))),xDep));
    
    % removes the day index dependency
    for i = 1:length(ii)
        pData.oP.yVar(ii(i)).xDep = xDep{ii(i)}(~strcmp(xDep{ii(i)},'iStimS'));
    end
    
    % determines which output variables have a day index dependency
    Stats = field2cell(pData.oP.yVar,'Stats');
    jj = find(cellfun(@(x)(any(strcmp(x,'iStimS'))),Stats));
    
    % removes the day index dependency
    for i = 1:length(jj)
        pData.oP.yVar(jj(i)).Stats = Stats{jj(i)}(~strcmp(Stats{jj(i)},'iStimS'));
    end    
end

% resets the group separation flag
pData.oP.sepGrp = false;

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
function p = resetPlotDataStruct(p0,grpStim,iApp)

% sets the fields not to be altered
pStr = {'T','iStim','gType'};

% retrieves the plotting data struct field names (removes the fields above)
fName = fieldnames(p0);
fName = fName(cellfun(@(x)(~any(strcmp(x,pStr))),fName));

% resets the static fields
p = struct();
for j = 1:length(pStr)
    eval(sprintf('p.%s = p0(1).%s;',pStr{j},pStr{j}));
end

% sets the data for the other fields
for j = 1:length(fName)
    % retrieves the values for the current field
    if (grpStim)
        Ynw = eval(sprintf('p0(iApp).%s{1};',fName{j}));
    else
        Ynw = cell2mat(cell2cell(field2cell(p0,fName{j}),0));
    end
    
    % sets the data values into the final plotting struct
    eval(sprintf('p.%s = Ynw;',fName{j}));    
end