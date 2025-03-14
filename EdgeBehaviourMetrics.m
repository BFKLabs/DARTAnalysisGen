% --- calculates the fly edge behaviour metrics over the duration
%     of an experiment (2D experiment only)
function pData = EdgeBehaviourMetrics(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Edge Behaviour (Metrics)';
pData.Type = {'Pop','Multi'};
pData.fType = [2 1 1 3];
pData.rI = initFuncReqInfo(pData);
pData.dcFunc = @dataCursorFunc;

% initialises the other fields  (if input argument provided)
if nargin == 1 
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
    [pData.hasSP,pData.hasRC] = deal(true);
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
rI.Shape = '2D (Circle)';
rI.Stim = 'None';
rI.Spec = 'None';
        
% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 11;                       
cP = setParaFields(nPara);

% sets the parameter list for list parameters
pList = cellfun(@num2str,num2cell([5 10 15 20 30]),'un',0);

% parameter tab strings
a = {'1 - General','2 - Thresholds'};

% sets the parameter fields
cP(1) = setParaFields(a{2},'Number',2,'vTol','Movement Threshold Speed (mm/s)',[0.1 3 false]);
cP(2) = setParaFields(a{2},'Number',3,'cTol','Circumferential Distance Threshold (mm)',[1 50 false]);
cP(3) = setParaFields(a{2},'Number',1,'tBefore','Pre-Edge Contact Duration (s)',[0.5 5 false]);
cP(4) = setParaFields(a{2},'Number',3,'tAfter','Post-Edge Contact Duration (s)',[1 10 false]);
cP(5) = setParaFields(a{2},'Number',0.5,'tMin','Minimum Edge Contact Duration (s)',[0 5 false]);
cP(6) = setParaFields(a{2},'Number',3,'dR','Outside Edge Distance (mm)',[0.1 5 false]);
cP(7) = setParaFields(a{1},'List',{2,pList},'dPhi','Approach Angle Bin Size (degrees)');
cP(8) = setParaFields(a{1},'Number',5,'nSm','Smoothing window size (frames)',[3 15 true]);
cP(9) = setParaFields(a{1},'Boolean',0,'useAll','Analyse Entire Experiment');
cP(10) = setParaFields(a{1},'Number',0,'T0','Start Time (min)',[0 inf true],{9,1});
cP(11) = setParaFields(a{1},'Number',10,'Tdur','Analysis Duration (min)',[1 inf true],{9,1});

% sets the tool-tip strings
cP(1).TTstr = 'The speed at which flies are considered to be moving';
cP(2).TTstr = 'The circumferential threshold distance for determining movement after edge contact';
cP(3).TTstr = 'The duration before contact with the wall used to determine approach angle';
cP(4).TTstr = 'The duration after contact with the wall used to determine edge behaviour';
cP(5).TTstr = 'Minimum duration fly has to be in the outer ring to be considered as edge contact';
cP(6).TTstr = 'The distance from outside edge whereby fly is considered in contact with wall';
cP(7).TTstr = 'Smoothing window size for the positional coordinates';
cP(8).TTstr = 'Angle bin size for the Approach Angle (Histogram) calculations';
cP(9).TTstr = 'Determines if all or part of the experiment will be used for the analysis';
cP(10).TTstr = 'Sets the start time of the analysis (minutes from the beginning)';
cP(11).TTstr = 'Duration of the sub-experiment analysis (in minutes)';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 8;
pP = setParaFields(nPara);

% sets the plot-type list
pList = {'Approach Angle (Histogram)','Approach Angle (Polar Plot)',...
         'Approach Angle (Vector)','Approach Angle (CDF)'...
         'Turn Direction (Histogram)','Turn Direction (Parameters)',...
         'Post-Contact (Movement Type)','Threshold Movement (Duration)',...
         'Post-Contact (Displacement)','Post-Contact (Duration)',...
         'Post-Contact (Direction Change)','Overall Fly Location'};
pList2 = {'Bar Graph','Boxplot'};  

% parameter tab strings
a = {'1 - General','2 - Histogram'};

% sets the parameter fields
pP(1) = setParaFields(a{1},'List',{1,pList},'pMet','Plot Metric');
pP(2) = setParaFields(a{1},'List',{1,pList2},'pType','Graph Plot Type',[],{1,[8:10 12]});
pP(3) = setParaFields(a{1},'Boolean',0,'plotGrid','Show Axis Gridlines',[],{1,[1 3:12]});
pP(4) = setParaFields(a{1},'Boolean',0,'plotAvg','Plot Average Value For Each Type',[],{1,7:11});
pP(5) = setParaFields(a{1},'Boolean',1,'plotErr','Show Error Bars/Outliers',[],{1,[3 6 8:10 12]});
pP(6) = setParaFields(a{2},'Boolean',0,'plotFit','Plot Fitted Curve',[],{1,[1 5]});
pP(7) = setParaFields(a{2},'Boolean',0,'showPara','Display Fitted Value',[],{1,5});
pP(8) = setParaFields(a{2},'Boolean',1,'usePDF','Plot Histogram As Proportional Values',[],{1,1});

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot,nApp)

% memory allocation
if nargin == 1; nApp = length(snTot.iMov.ok); end
pF = setFormatFields(nApp);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],'',1);
pF.yLabel = setFormatFields([],'',1);
pF.Legend = setFormatFields(setupFontStruct('FontSize',12),'',1);
pF.Axis = setFormatFields([],[]);

% sets the apparatus names as the titles
if nargin == 1
    for i = 1:nApp
        pF.Title(i).String = snTot.iMov.pInfo.gName{i};
    end
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot,false);
[Stats1,Stats2,Stats3] = deal({'TTest'},{'CompMulti','phiB'},{'Comp'});
[Type1,Type2,Type3] = deal(4,3,2);
[xDep1,xDep2,xDep3] = deal({'phiB'},{'xCDF'},{'phiB','scmStr'});

% sets the independent output variables
oP = addXVarField(oP,'Approach Angle Bin','phiB','Group');
oP = addXVarField(oP,'Activity Type','scmStr','Group');
oP = addXVarField(oP,'Approach CDF Angles','xCDF','Other');

% sets the dependent output variables
oP = addYVarField(oP,'Approach Histogram (Raw)','PA',[],Type1,xDep1);
oP = addYVarField(oP,'Approach Histogram (Fit)','PAfit',[],Type1,xDep1);
oP = addYVarField(oP,'Turn Direction (Raw)','TD',[],Type1,xDep1);
oP = addYVarField(oP,'Turn Direction (Fit)','TDfit',[],Type1,xDep1);
oP = addYVarField(oP,'Approach CDF Proportions','fCDF',[],Type1,xDep2);
oP = addYVarField(oP,'Slope Factor','nTD',Stats1,[],[]);
oP = addYVarField(oP,'Slope Factor (Mean)','nTD_mn',[],Type3,[]);
oP = addYVarField(oP,'Slope Factor (SEM)','nTD_sem',[],Type3,[]);
oP = addYVarField(oP,'Shape Factor','ATD',Stats1,[],[]);
oP = addYVarField(oP,'Shape Factor (Mean)','ATD_mn',[],Type3,[]);
oP = addYVarField(oP,'Shape Factor (SEM)','ATD_sem',[],Type3,[]);
oP = addYVarField(oP,'Movement Type','DM',[],Type1,xDep3);
oP = addYVarField(oP,'Outer Region %age','PrE',Stats3,Type2,[],1);
oP = addYVarField(oP,'Outer Region %age (Mean)','PrE_mn',[],Type3,[]);
oP = addYVarField(oP,'Outer Region %age (SEM)','PrE_sem',[],Type3,[]);
oP = addYVarField(oP,'Outer Region %age (Median)','PrE_md',[],Type3,[]);
oP = addYVarField(oP,'Movement Duration','TM',Stats2,[],xDep1);
oP = addYVarField(oP,'Movement Duration (Mean)','TM_mn',[],Type1,xDep1);
oP = addYVarField(oP,'Movement Duration (SEM)','TM_sem',[],Type1,xDep1);
oP = addYVarField(oP,'Movement Duration (Median)','TM_md',[],Type1,xDep1);
oP = addYVarField(oP,'Edge Duration','TE',Stats2,[],xDep1);
oP = addYVarField(oP,'Edge Duration (Mean)','TE_mn',[],Type1,xDep1);
oP = addYVarField(oP,'Edge Duration (SEM)','TE_sem',[],Type1,xDep1);
oP = addYVarField(oP,'Edge Duration (Median)','TE_md',[],Type1,xDep1);
oP = addYVarField(oP,'Edge Displacement','DE',Stats2,[],xDep1);
oP = addYVarField(oP,'Edge Displacement (Mean)','DE_mn',[],Type1,xDep1);
oP = addYVarField(oP,'Edge Displacement (SEM)','DE_sem',[],Type1,xDep1);
oP = addYVarField(oP,'Edge Displacement (Median)','DE_md',[],Type1,xDep1);

% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
pP = retParaStruct(dcObj.pData.pP);
sP = retParaStruct(dcObj.pData.sP);

% other initialisations
phiB = dcObj.plotD{1}(1).phiB;
grpName = dcObj.pData.appName(sP.Sub.isPlot);

% sets the common class fields
dcObj.yName = pP.pMet; 
dcObj.grpName = grpName;
dcObj.useGrpHdr = true;
dcObj.useXGrp = false;
dcObj.combFig = false;

% sets the metric specific fields
switch pP.pMet
    case 'Approach Angle (Polar Plot)'                
        % case is the polar plot metrics        
        PA = field2cell(dcObj.plotD{1}(sP.Sub.isPlot),'PA');
        iAx = dcObj.getSelectAxesIndex;
        
        % sets the class object fields
        dcObj.pType = 'Polar';
        dcObj.xName = 'Approach Angle';    
        dcObj.yName = 'Approach Angle (Counts)'; 
        [dcObj.xGrp,dcObj.yGrp] = deal(phiB,PA{iAx});
    
    case 'Approach Angle (Vector)'
        % case is the boxplot only metrics
        
        % sets the class object fields
        dcObj.pType = 'Boxplot';
        [dcObj.xGrp,dcObj.useGrpHdr] = deal(grpName,false);
        [dcObj.xName,dcObj.yUnits] = deal('Group Name','unitless');         
        
    case {'Threshold Movement (Duration)','Post-Contact (Displacement)',...
          'Post-Contact (Duration)','Overall Fly Location'}
        % case is the bar/boxplot metrics
        
        % sets the class object fields
        dcObj.pType = pP.pType;
        dcObj.xName = 'Approach Angle';
        dcObj.xGrp = phiB;
        
        % sets the y-axis units
        switch pP.pMet
            case 'Post-Contact (Displacement)'
                % case is the post-contact displacement
                dcObj.yUnits = 'mm';
                
            case {'Threshold Movement (Duration)',...
                  'Post-Contact (Duration)'}
                % case is the threshold movement (duration)
                dcObj.yUnits = 'sec';
                
            case 'Overall Fly Location'
                % case is the duration specific metrics
                dcObj.yUnits = '%';
                dcObj.xGrp = grpName;
                dcObj.xName = 'Group Name';
                dcObj.useGrpHdr = false;
        end
        
    case 'Post-Contact (Movement Type)'
        % case is the stacked barplot metrics        
        
        % sets the class object fields
        dcObj.yUnits = '%';
        dcObj.pType = 'Stacked Bar Graph';          
        dcObj.yUnits = 'unitless';
        dcObj.xName = 'Approach Angle';        
        dcObj.xGrp2 = {'Stationary','Combination','Moving'};
        
        if pP.plotAvg
            dcObj.xGrp = grpName;
            dcObj.useGrpHdr = false;
            dcObj.combFig = true;
        else
            dcObj.xGrp = phiB;
        end
        
    case 'Approach Angle (CDF)'
        % case is the trace only metrics
        
        % sets the class object fields
        dcObj.yUnits = '%';
        dcObj.pType = 'Trace';
        dcObj.xName = 'Approach Angle';
        dcObj.xGrp = phiB;
        dcObj.yUnits = 'unitless';
        dcObj.combFig = true;
        
    case {'Approach Angle (Histogram)','Turn Direction (Duration)'}
        % case is the bar + trace metrics

        % sets the common class fields
        dcObj.xName = 'Approach Angle';
        dcObj.xGrp = phiB;
        
        % sets the class object fields
        if strcmp(get(dcObj.evnt.Target,'Type'),'bar')
            % sets the common field types
            dcObj.pType = 'Bar Graph';
            dcObj.yUnits = 'Frequency';
            
        else
            % sets the common field types
            dcObj.pType = 'Trace';
            dcObj.useXGrp = true;
            dcObj.tUnits = 'degrees';
            dcObj.yUnits = 'unitless';
        end
        
    otherwise
        % case is the bar graph only metrics
        dcObj.pType = 'Bar Graph';
        
        % sets the class object fields
        switch pP.pMet
            case 'Turn Direction (Parameters)'
                % case is the turn direction parameters
                dcObj.grpName = {'Slope Factor','Shape Factor',...
                                 'Correlation'};
                dcObj.xGrp = grpName;
                dcObj.xName = 'Approach Angle';  
                
            case 'Post-Contact (Direction Change)'
                % case is the post-contact direction change
                dcObj.yUnits = '%';                
                dcObj.xGrp = phiB;
        end
end

% %

% pList2 = {'Bar Graph','Boxplot'}; 

% sets up the data cursor string
dTxt = dcObj.setupCursorString();

% ----------------------------------------------------------------------- %
% ---                       CALCULATION FUNCTION                      --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [plotD,ok] = calcFunc(snTot,pData,gPara,cP,varargin)

% initialises the calculation parameters (if not already initialised)
if nargin == 3
    % retrieves the parameter struct
    cP = retParaStruct(pData.cP,gPara);
end

% checks to see if the solution struct has the sub-region data struct
snTotL = snTot(1);
ok = checkFuncPara({'HasSubRegionStruct'},cP,snTotL);
if ~ok; plotD = []; return; end

% sets the movement calculation type
cP.movType = 'Absolute Speed';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);
scmStr = {'Stationary','Combined','Moving'};

% calculation parameters
dPhi = str2double(cP.dPhi);
phiT = (-(90+dPhi):dPhi:(90+dPhi));
phiB = phiT(2:(end-1));

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                         'PA',[],'DT',[],'TD',[],'DM',[],...
                         'phiB',phiB,'iBin',[],'PAfit',[],...
                         'TDfit',[],'pTD',[],'ATD',[],'nTD',[],'R2TD',[],...                                 
                         'TM',[],'TM_mn',[],'TM_sem',[],'TM_md',[],...
                         'TE',[],'TE_mn',[],'TE_sem',[],'TE_md',[],...
                         'DE',[],'DE_mn',[],'DE_sem',[],'DE_md',[],...                                 
                         'PrE',[],'PrE_mn',[],'PrE_sem',[],'PrE_md',[],...
                         'xCDF',[],'fCDF',[],'scmStr',scmStr);                        

% memory allocation
[PhiApp,Dmove,Tdir,Dturn,Tmove,Dedge,Tedge] = deal(cell(1,nApp));

% creates the waitbar figure
wStr = {'Turn Detection','Turning Metric Calculations'};
h = ProgBar(wStr,'Turning Behaviour Calculations');

% ------------------------------------------------------- %
% --- INTER-FRAME DISTANCE CALCULATION & THRESHOLDING --- %
% ------------------------------------------------------- %
    
% sets the time arrays for each experiment
T = cellfun(@(x)(cell2mat(x)),field2cell(snTot,'T'),'un',0);

% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure 
    wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
    if h.Update(1,wStrNw,i/(1+nExp))
        [plotD,ok] = deal([],false);
        return
    else
        h.Update(2,wStr{2},0);
    end

    % calculates the video frame rate and experiment apparatus indices        
    iApp = find(~cellfun('isempty',snTot(i).iMov.flyok));     

    % sets the relevant time points and apparatus indices for this expt
    if cP.useAll
        % uses all the time points
        ii = 1:length(T{i});
        nFrm = length(ii);
    else
        % use only the points from the start to the duration end
        ii = (T{i} >= 60*cP.T0) & (T{i} <= 60*(cP.T0 + cP.Tdur));
        nFrm = sum(ii);
    end    
    
    % determines the number of search frames for the pre/post edge contact
    % index search
    dTmn = median(diff(T{i}));
    [nFrmB,nFrmA] = deal(ceil(cP.tBefore/dTmn),ceil(cP.tAfter/dTmn));
    Dtol = cP.vTol*dTmn;
    
    % sets the relevant x/y-locations for the current experiment 
    [dPx,dPy,R] = get2DCoordsBG(snTot(i),iApp,ii);
    if length(iApp) == 1
        [dPx,dPy,R] = deal({dPx},{dPy},{R});
    end
        
    % determines the turning events for each frame
    for j = 1:length(iApp)
        % updates the waitbar figure 
        wStrNw = sprintf('%s (Apparatus %i of %i)',wStr{2},j,length(iApp));
        if h.Update(2,wStrNw,j/length(iApp))
            [plotD,ok] = deal([],false);
            return
        end        
        
        % calculates the angle of the fly wrt the circle centre
        PhiC = NaN(size(dPx{j}));
        for k = 1:size(dPx{j},2)            
            % smooths the non-NaN components of the angle array
            kGrp = getGroupIndex(~isnan(dPx{j}(:,k)));
            for kk = 1:length(kGrp)
                dPx{j}(kGrp{kk},k) = smooth(dPx{j}(kGrp{kk},k),cP.nSm);
                dPy{j}(kGrp{kk},k) = smooth(dPy{j}(kGrp{kk},k),cP.nSm);                                
            end
            
            % angle calculation
            PhiC(:,k) = unwrap(atan2(dPy{j}(:,k),dPx{j}(:,k)));            
        end          
        
        % calculates the distance from the circle centre
        sFac = snTot(i).sgP.sFac;
        [Rho,dR] = deal(sqrt(dPx{j}.^2 + dPy{j}.^2),R{j}-cP.dR/sFac);
                 
        % determines the time points where the flies have moved into the
        % edge region. from this calculate the proportion of time that the
        % flies occupy the outer region
        onEdge = Rho > repmat(dR(:)',nFrm,1);        
        [onEdge,a] = deal(num2cell(onEdge,1),cell(length(dR),1));
                
        % determines the index groups where the fly is in contact with the
        % edge of the arena. groupings 
        iGrp = cellfun(@(x)(getGroupIndex(x)),onEdge,'un',0);                        
        [PhiAppNw,DmoveNw,TmoveNw,TdirNw,DedgeNw,TedgeNw,DturnNw] = deal(a);
        for k = 1:length(iGrp)
            % memory allocations
            if ~isempty(iGrp{k})
                % determines the start time for contact with the wall, and
                % the duration that the fly spent in the outer edge region
                T0 = cellfun(@(x)(x(1)),iGrp{k});
                Tf = cellfun(@(x)(x(end)),iGrp{k});
                
                % calculates the within and between group times
                iWG = cellfun(@(x)(diff(T{i}(x([1 end])))),iGrp{k});
                iBG = [inf;(T0(2:end)-Tf(1:end-1))];
                
                % removes the groups that either within the frame tolerance
                % of the previous group or are too short. also removes
                % groups that have a start time less than tBefore, or is
                % less than tAfter from the end of the time vector
                ii = (iWG > ceil(cP.tMin/dTmn)) & (iBG > nFrmB) & ...
                     (T0 > nFrmB) & (T0 < (nFrm - nFrmA));
                
                % removes any groups that have start indices within the
                % frame tolerance of the preceding group
                if any(ii)
                    % ----------------------------- %                    
                    % --- INVALID GROUP REMOVAL --- %
                    % ----------------------------- %                    
                    
                    % removes all the groups that have edge contact times
                    % that are too close to the previous event
                    iPr = find(ii,1,'first');                                
                    for iNw = (iPr+1):length(iGrp{k})
                        % only check valid groups
                        if ii(iNw)
                            % determines if the initial index of the new group
                            % relative to the previous is greater than the
                            % frame tolerance
                            if (T0(iNw) - T0(iPr)) > (nFrmB+nFrmA)
                                % if so, then update the previous frame index
                                % to the new index
                                iPr = iNw;
                            else
                                % otherwise, flag the new group to be invalid
                                ii(iNw) = false;
                            end
                        end
                    end

                    % removes any of the non-valid groups
                    iGrp{k} = iGrp{k}(ii);                    
                    
                    % sets the pre/post edge contact index arrays
                    iB = cellfun(@(x)(x(1)+(-nFrmB:0)),iGrp{k},'un',0);
                    iA = cellfun(@(x)(x(1)+(0:nFrmA)),iGrp{k},'un',0);

                    % ---------------------------------- %                    
                    % --- APPROACH ANGLE CALCULATION --- %
                    % ---------------------------------- %
                    
                    % calculates the approach angle for all events where
                    % the fly enters the outside region
                    PhiAppNw{k} = cellfun(@(x)(calcApproachAngle(...
                                    dPx{j}(x,k),dPy{j}(x,k))*180/pi),iB);       
                                
                    % ---------------------------------------------- %                    
                    % --- CIRCUMFERENTIAL DISTANCE MOVEMENT TIME --- %
                    % ---------------------------------------------- %                    
                                
                    % sets the angles
                    PhiA = cell2mat(cellfun(@(x)(PhiC(x(1),k) - ...
                                    PhiC(x,k)'),iA,'un',0));
                                
                    % determines the indices where the fly moves
                    % circumferentially through a distance, cTol
                    cTol = cP.cTol/dR(k);
                    A = sign(PhiA).*(abs(PhiA) > cTol);
                    iMove = cellfun(@(x)(find([x,1],1,'first')-1),...
                                    num2cell(A,2));

                    % calculates the change in the circumferential angle
                    % over the duration of the sub-sequence.
                    dPhiA = diff(PhiA,[],2);
                    dPhiA = dPhiA.*(abs(dPhiA) > Dtol/dR(k));                                
                                
                    % calculates the movement time (removes any flies that
                    % did not make it to the threshold distance)
                    TmoveNw{k} = iMove*dTmn;
                    TmoveNw{k}(iMove > nFrmA) = NaN;
                    
                    % sets the turn direction for each fly                     
                    ind = sub2ind(size(A),(1:size(A,1))',min(iMove+1,nFrmA+1));
                    TdirNw{k} = A(ind);
                    
                    % calculates if fly performed any double turns
                    DturnNw{k} = any(sign(dPhiA) == -1,2) & any(sign(dPhiA) == 1,2);    
                    
                    % ------------------------------------- %                    
                    % --- EDGE BEHAVIOUR CLASSIFICATION --- %
                    % ------------------------------------- %                    
                    
                    % calculates the inter-frame displacement of the flies
                    % during the after edge contact phase. from this,
                    % determine the average frame movement of the flies
                    % during this phase (average of the boolean flags).
                    % removes any flies that did not make the
                    % circumferential distance threhold (sets value to 0)
                    DedgeT = cellfun(@(x)(sqrt(diff(dPx{j}(x,k)).^2 + ...
                            diff(dPy{j}(x,k)).^2)'),iA,'un',0);                           
                    DedgeNw{k} = sum(cell2mat(DedgeT),2);
                    DmoveNw{k} = mean(cell2mat(DedgeT) > Dtol,2);
                    DmoveNw{k}(isnan(TmoveNw{k})) = 0;
                    
                    % determines if duration at which the fly stays within
                    % the edge region
                    TedgeNw{k} = cellfun('length',iGrp{k});
                end            
            end
        end       
        
        % append that data to the overall arrays
        ii = ~cellfun('isempty',PhiAppNw);
        PhiApp{iApp(j)} = [PhiApp{iApp(j)};cell2mat(PhiAppNw(ii))];
        Dmove{iApp(j)} = [Dmove{iApp(j)};cell2mat(DmoveNw(ii))];        
        Dturn{iApp(j)} = [Dturn{iApp(j)};cell2mat(DturnNw(ii))];        
        Tdir{iApp(j)} = [Tdir{iApp(j)};cell2mat(TdirNw(ii))];        
        Dedge{iApp(j)} = [Dedge{iApp(j)};cell2mat(DedgeNw(ii))];
        Tedge{iApp(j)} = [Tedge{iApp(j)};cell2mat(TedgeNw(ii))*dTmn];                
        Tmove{iApp(j)} = [Tmove{iApp(j)};cell2mat(TmoveNw(ii))];
        
        % sets the 
        [kk,indF] = deal(iApp(j),1:length(onEdge));
        plotD(kk).PrE(1,indF,i) = cellfun(@(x)(100*mean(x)),onEdge,'un',0);
%         [TM{i,kk},DE{i,kk}] = deal(TmoveNw,DedgeNw);
%         TE{i,kk} = cellfun(@(x)(x*dTmn),TedgeNw,'un',0);  
    end   
end    

% sets/calculates the raw, mean and SEM proportional activity values
for i = 1:nApp   
    % updates the waitbar figure 
    wStrNw = sprintf('Setting Final Values (Apparatus %i of %i)',i,nApp);
    if h.Update(2,wStrNw,i/nApp)
        [plotD,ok] = deal([],false);
        return
    end     
    
    % removes all the NaN values from the arrays
    ii = ~isnan(PhiApp{i});
    [PhiApp{i},Dturn{i}] = deal(PhiApp{i}(ii),Dturn{i}(ii)*100);    
    Tdir{i} = Tdir{i}(ii);
                
    % sets the binned angles array    
    phiTB = phiT + (dPhi/2); phiTB(end) = inf; 
    iBin = cellfun(@(x)(find(x<phiTB,1,'first')),num2cell(PhiApp{i}));            
    plotD(i).iBin = cellfun(@(x)(find(iBin==x)),...
                    num2cell(2:length(phiTB)-1),'un',0);
                        
    % fits the parameters for the approach angle PDF 
    PA = cellfun(@(x)(PhiApp{i}(x)),plotD(i).iBin,'un',0);    
    plotD(i).PA = cellfun('length',PA);
    [plotD(i).pPA,plotD(i).PAfit] = fitBimodalDist(plotD(i).phiB,plotD(i).PA);                          
                        
    % fits the parameters for the turning direction distributions
    TD = cellfun(@(x)(Tdir{i}(x)),plotD(i).iBin,'un',0);
    plotD(i).TD = cellfun(@(x)(mean(x,'omitnan')),TD);    
    [plotD(i).pTD,plotD(i).TDfit] = ...
                fitTurnDirDist(plotD(i).phiB,plotD(i).TD,plotD(i).PA);                                        

    % sets the fitted Turn Direction (Parameters)
    [plotD(i).nTD,plotD(i).ATD,plotD(i).R2TD] = ...
                    deal(plotD(i).pTD.n,plotD(i).pTD.A,plotD(i).pTD.R2);
                
    % sets the other metrics into the plot data array    
    plotD(i).TM = cellfun(@(x)(Tmove{i}(x)'),plotD(i).iBin,'un',0);       
    plotD(i).DT = cellfun(@(x)(Dturn{i}(x)'),plotD(i).iBin,'un',0);       
    plotD(i).DE = cellfun(@(x)(Dedge{i}(x)'),plotD(i).iBin,'un',0);       
    plotD(i).TE = cellfun(@(x)(Tedge{i}(x)'),plotD(i).iBin,'un',0);       
    
    % calculates the movement type proportions (moving, stationary or both)
    DM = cellfun(@(x)(Dmove{i}(x)),plotD(i).iBin,'un',0);
    plotD(i).DM = zeros(length(DM),3);
    plotD(i).DM(:,1) = 100*cellfun(@(x)(sum(x==0)/length(x)),DM);
    plotD(i).DM(:,3) = 100*cellfun(@(x)(sum(x==1)/length(x)),DM);
    plotD(i).DM(:,2) = 100 - sum(plotD(i).DM(:,[1 3]),2);
             
    % calculates the metric statistics
    plotD(i) = calcMetricStats(plotD(i),'PrE');     
    plotD(i) = calcMetricStats(plotD(i),'DE',6);  
    plotD(i) = calcMetricStats(plotD(i),'TE',6);
    plotD(i) = calcMetricStats(plotD(i),'TM',6);  
    
    % calculates the Approach Angle (CDF) distributions
    [plotD(i).xCDF,plotD(i).fCDF] = calcAngleCDF(PhiApp{i});    
end
  
% closes the waitbar figure
if ~h.Update(1,'Activity Calculations Complete!',1)
    if nargin < 5; h.closeProgBar(); end
end

% ----------------------------------------------------------------------- %
% ---                        PLOTTING FUNCTION                        --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,varargout] = plotFunc(snTot,pData,plotD,ind)

% global variables
global regSz newSz

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
nApp = length(ind); if nApp == 0; return; end
p = plotD{1}(ind);

% sets the flag determining if there are going to be subplots
[isSub,setYLim,rotXAx,setXTick] = deal(true,true,false,true);
[dPhi,pStep] = deal(str2double(cP.dPhi),30);

% sets the x-axis tick marks and strings
xTickLbl = cellfun(@num2str,num2cell(-90:pStep:90),'un',0);
xTick = (0:(pStep/dPhi):(180/dPhi)) + 1;

% other parameters
[pDel,lgOfs] = deal(0.1,0.00);

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %
   
% initialisations
lgHorz = false;
[Yfit,Ysem,lStr,fitVal,fitLoc,yLim] = deal([]);

% sets the default x-axis plot values, limits and axis strings
X = 1:length(p(1).phiB);
xLim = X([1 end]) + 0.5*[-1 1];
[xStr,pType] = deal('Approach Angle (deg)','bar');

% sets the plot values/axis labels based on the metric type
switch pP.pMet
    case ('Approach Angle (Histogram)')
        % sets the plot values    
        Y = field2cell(p,'PA');
        N = cellfun(@(x)(sum(x)),Y,'un',0);
        if pP.usePDF
            % if using the the PDF values, then scale them             
            Y = cellfun(@(x,y)(x/y),Y,N,'un',0); 
        end
        
        % sets the fitted values (if required)
        if pP.plotFit
            Yfit = field2cell(p,'PAfit');            
            if ~pP.usePDF
                % if not using the the PDF values, then scale them
                Yfit = cellfun(@(x,y)(x*y),Yfit,N,'un',0); 
            end
        end
            
        % sets the x/y-axis label strings        
        if pP.usePDF; yStr = 'Proportion'; else yStr = 'Frequency'; end
        
    case ('Approach Angle (Polar Plot)')
        % sets the plot data
        [phi,pType,yStr] = deal(p(1).phiB,'polar','');
        [dphi,d2r,pF.Axis.Font.Color] = deal(diff(phi([1 2])),(pi/180),'w');
        
        % sets the plot data
        Y = cellfun(@(x)(x/max(x)),field2cell(p,'PA'),'un',0);               
        
    case ('Approach Angle (Vector)')
        %
        [rotXAx,isSub,xStr] = deal(true,false,'');
        [pType,yStr] = deal('boxplot','Vector Length (Normalised)');
        Y = combineNumericCells(cellfun(@(x)(x/max(x)),field2cell(p,'PA'),'un',0));        
        
        % sets the x-axis properties        
        [X,xTick] = deal(1:length(p));
        xTickLbl = field2cell(pF.Title(ind),'String');
        xLim = xTick([1 end]) + 0.5*[-1 1];        
        
        % sets the title/x-label properties
        pF.Title(1).String = 'Approach Angle Vector Lengths';
        pF.xLabel(1).String = '';
        pF.xLabel(1).Font.FontSize = 1;
        
    case ('Approach Angle (CDF)')
        % sets the plot values  
        [isSub,setXTick,pType] = deal(false,false,'plot');
        [X,Y] = deal(p(1).xCDF,100*cell2mat(field2cell(p','fCDF')));
        [xLim,xTick] = deal([-90 90],-90:10:90);
        
        % set main title and y-axis strings
        pF.Title(1).String = 'Approach Angle (CDF)';
        yStr = 'Cumulative Percentage'; 
        lStr = snTot(1).iMov.pInfo.gName(ind);
        
%     case ('Approach Angle Parameters')
%         % sets the plot values  
%         [X,nApp,m,n,setYLim,rotXAx] = deal(1:nApp,4,1,4,false,true);
%         [Y,Ysem] = deal(cell(nApp,1));
%         [pVal,ind] = deal(cell2mat(field2cell(p,'pPA')),1:nApp);
%         
%         % sets the x-axis properties
%         [xTick,xTickLbl] = deal(1:length(p),field2cell(pF.Title,'String'));
%         xLim = xTick([1 end]) + 0.5*[-1 1];
%         pF = initPlotFormat(snTot(1),nApp);
%         
%         % resets the title strings
%         [tStr,fStr] = deal({'\mu','\sigma','A','R^2'},{'mu','sigma','A','R2'});        
%         for i = 1:nApp
%             % sets the sub-plot title strings
%             pF.Title(i).String = tStr{i};
%             [pF.xLabel(i).String,pF.yLabel(i).String] = deal('');
%             if (i > 1)
%                 pF.xLabel(i).Font = pF.xLabel(1).Font;
%                 pF.yLabel(i).Font = pF.yLabel(1).Font;
%             end
%                         
%             % sets the plot values and SEM values (if they exist)
%             Ynw = cell2mat(field2cell(pVal,fStr{i}));
%             Y{i} = Ynw(:,1);
%             if (size(Ynw,2) == 2); Ysem{i} = Ynw(:,2); end
%         end
%                 
%         % sets the y-axis label strings        
%         [yStr,xStr] = deal('Parameter Value','Group Name');
        
    case ('Turn Direction (Histogram)')
        % sets the plot values    
        [fitVal,fitLoc] = deal(field2cell(p','pTD'),'NorthEast');
        Y = field2cell(p,'TD')';
        
        % sets the fitted values (if required)
        if pP.plotFit
            Yfit = field2cell(p,'TDfit');            
        end
            
        % sets the y-axis label strings        
        yStr = 'Turn Direction';
        
    case ('Turn Direction (Parameters)')
        % sets the plot values  
        [X,nApp,m,n,setYLim,rotXAx,ii] = deal(1:nApp,3,1,3,false,true,ind);
        [Y,Ysem] = deal(cell(nApp,1));
        ind = 1:nApp;        
                        
        % sets the x-axis properties
        [xTick,xTickLbl] = deal(1:length(p),field2cell(pF.Title(ii),'String'));
        xLim = xTick([1 end]) + 0.5*[-1 1];
        pF = initPlotFormat(snTot(1),nApp);
        
        % resets the title strings
        [tStr,fStr] = deal({'n','A','R^2'},{'nTD','ATD','R2TD'});
        for i = 1:nApp
            % sets the sub-plot title strings
            pF.Title(i).String = tStr{i};
                        
            % sets the plot values and SEM values (if they exist)
            Ynw = cell2mat(field2cell(p,fStr{i}));
            Y{i} = Ynw(:,1);
            if size(Ynw,2) == 2; Ysem{i} = Ynw(:,2); end
        end
                
        % sets the y-axis label strings        
        [yStr,xStr] = deal('Parameter Value','');
                   
    case ('Post-Contact (Movement Type)')
        % retrieves 
        pType = 'barstack';
        Y = field2cell(p,'DM');       
        if pP.plotAvg
            % combines the data into a single array
            [rotXAx,isSub] = deal(true,false);
            Y = cell2mat(cellfun(@(x)(mean(x,1,'omitnan')),Y,'un',0));
                        
            % sets the x-axis properties
            [X,xTick] = deal(1:length(p));
            [xTickLbl,xStr] = deal(field2cell(pF.Title(ind),'String'),'');
            xLim = xTick([1 end]) + 0.5*[-1 1];            
            
            % resets the title string
            pF.Title(1).String = 'Post-Contact (Movement Type)';
        end
        
        % set legend strings and the y-axis string/limits
        [lgLoc,lgHorz] = deal('EastOutside',true);
        lStr = {'Stationary','Combination','Moving'}';
        [yStr,yLim] = deal('Percentage',[0 100]);        
                
    case ('Threshold Movement (Duration)')
        % sets the plot values
        TM = field2cell(p,'TM');
        if strcmp(pP.pType,'Boxplot')
            % sets the plot type
            pType = 'boxplot';            
            if pP.plotAvg
                % combines the data into a single array for each group
                TM = cellfun(@(x)(cell2mat(x(:)')),TM,'un',0);                                 
                [Y,isSub,rotXAx] = deal(combineNumericCells(TM'),false,true);
            else
                % combines the data into a single array for each time group
                Y = cellfun(@(x)(combineNumericCells(x')),TM,'un',0);
            end
        else
            % sets the plot type
            pType = 'bar';
            if pP.plotAvg
                % combines the data into a single cell array
                TM = cellfun(@(x)(cell2mat(x(:)')),TM,'un',0);                
                
                % recalculates the mean/SEM values                
                [isSub,rotXAx] = deal(false,true);
                Y = cellfun(@(x)(mean(x,'omitnan')),TM);
                Ysd = cellfun(@(x)(std(x,[],'omitnan')),TM);
                Ysem = Ysd./sqrt(cellfun('length',TM));
            else                
                [Y,Ysem] = field2cell(p,{'TM_mn','TM_sem'});
            end
        end
        
        % resets the titles/labels and axis properties
        if pP.plotAvg
            % sets the x-axis properties
            [X,xTick] = deal(1:length(p));
            [xTickLbl,xStr] = deal(field2cell(pF.Title(ind),'String'),'');
            xLim = xTick([1 end]) + 0.5*[-1 1];            
            
            % resets the title string
            pF.Title(1).String = 'Distance Movement Duration';            
        end
        
        % sets the y-axis string
        yStr = 'Threshold Time (sec)';                    
        
    case ('Post-Contact (Displacement)')
        % sets the plot values
        DE = field2cell(p,'DE');
        if strcmp(pP.pType,'Boxplot')
            % sets the plot type
            pType = 'boxplot';            
            if pP.plotAvg
                % combines the data into a single array for each group
                DE = cellfun(@(x)(cell2mat(x(:)')),DE,'un',0);                                 
                [Y,isSub,rotXAx] = deal(combineNumericCells(DE'),false,true);
            else
                % combines the data into a single array for each time group
                Y = cellfun(@(x)(combineNumericCells(x')),DE,'un',0);
            end
        else
            % sets the plot type
            pType = 'bar';
            if pP.plotAvg
                % combines the data into a single cell array
                DE = cellfun(@(x)(cell2mat(x(:)')),DE,'un',0);                
                
                % recalculates the mean/SEM values
                [isSub,rotXAx] = deal(false,true);
                Y = cellfun(@(x)(mean(x,'omitnan')),DE);
                Ysd = cellfun(@(x)(std(x,[],'omitnan')),DE);
                Ysem = Ysd./sqrt(cellfun('length',DE));
            else                
                [Y,Ysem] = field2cell(p,{'DE_mn','DE_sem'});
            end                                    
        end
        
        % resets the titles/labels and axis properties
        if pP.plotAvg
            % sets the x-axis properties
            [X,xTick] = deal(1:length(p));
            [xTickLbl,xStr] = deal(field2cell(pF.Title(ind),'String'),'');
            xLim = xTick([1 end]) + 0.5*[-1 1];            
            
            % resets the title string
            pF.Title(1).String = 'Edge Region Displacement';            
        end        
        
        % sets the y-axis string        
        yStr = 'Distance (mm)';
        
    case ('Post-Contact (Duration)')
        % sets the plot values
        TE = field2cell(p,'TE');
        if strcmp(pP.pType,'Boxplot')
            % sets the plot type
            pType = 'boxplot';            
            if pP.plotAvg
                % combines the data into a single array for each group
                TE = cellfun(@(x)(cell2mat(x(:)')),TE,'un',0);                                 
                [Y,isSub,rotXAx] = deal(combineNumericCells(TE'),false,true);
            else
                % combines the data into a single array for each time group
                Y = cellfun(@(x)(combineNumericCells(x')),TE,'un',0);
            end
        else
            % sets the plot type
            pType = 'bar';
            if pP.plotAvg
                % combines the data into a single cell array
                TE = cellfun(@(x)(cell2mat(x(:)')),TE,'un',0);                                 
                
                % recalculates the mean/SEM values
                [isSub,rotXAx] = deal(false,true);
                Y = cellfun(@(x)(mean(x,'omitnan')),TE);
                Ysd = cellfun(@(x)(std(x,[],'omitnan')),TE);
                Ysem = Ysd./sqrt(cellfun('length',TE));
            else                
                [Y,Ysem] = field2cell(p,{'TE_mn','TE_sem'});
            end                                    
        end
        
        % resets the titles/labels and axis properties
        if pP.plotAvg
            % sets the x-axis properties
            [X,xTick] = deal(1:length(p));
            [xTickLbl,xStr] = deal(field2cell(pF.Title(ind),'String'),'');
            xLim = xTick([1 end]) + 0.5*[-1 1];            
            
            % resets the title string
            pF.Title(1).String = 'Edge Region Duration';            
        end        
        
        % sets the y-axis string        
        yStr = 'Duration (s)';                       
        
    case ('Post-Contact (Direction Change)')
        % sets the plot values
        [DT,pType] = deal(field2cell(p,'DT'),'bar');
        if pP.plotAvg
            % combines the data into a single cell array
            DT = cellfun(@(x)(cell2mat(x(:)')),DT,'un',0);                

            % recalculates the mean/SEM values
            [isSub,rotXAx] = deal(false,true);
            Y = cellfun(@(x)(mean(x,'omitnan')),DT);
        else                
            Y = cellfun(@(x)(cellfun(@(y)(mean(y,'omitnan')),x)),DT,'un',0);
        end                                    
        
        % resets the titles/labels and axis properties
        if pP.plotAvg
            % sets the x-axis properties
            [X,xTick] = deal(1:length(p));
            [xTickLbl,xStr] = deal(field2cell(pF.Title(ind),'String'),'');
            xLim = xTick([1 end]) + 0.5*[-1 1];            
            
            % resets the title string
            pF.Title(1).String = 'Post-Contact (Direction Change)';            
        end        
        
        % sets the y-axis string        
        yStr = 'Direction Change (%)';        
        
    case ('Overall Fly Location')
        % sets the plot values
        [X,isSub,rotXAx] = deal(1:length(p),false,true);
        if strcmp(pP.pType,'Boxplot')
            [pType,yLim] = deal('boxplot',[0 100]+5*[-1 1]);
            Y = combineNumericCells(field2cell(p,'PrE')');
        else
            [pType,yLim] = deal('bar',[0 100]);
            [Y,Ysem] = field2cell(p,{'PrE_mn','PrE_sem'},1);
        end        
        
        % sets the x-axis tick indices/labels and limits
        [xTick,xTickLbl] = deal(X,snTot(1).iMov.pInfo.gName(ind));
        xLim = X([1 end]) + 0.5*[-1 1];
        
        % sets the y-axis string
        pF.Title(1).String = 'Fly Location In Outer Region';
        pF.xLabel(1).Font.FontSize = 1;
        [xStr,yStr] = deal('','Percentage of Time');
end

%
if isSub
    % reformats the font data struct for the current subplot count
    if isempty(m); szMx = 1; else; szMx = max([m n]); end
    pF = retFormatStruct(pF,szMx);
    
    % sets the y-label indices
    [hAx,hPlt] = deal(zeros(nApp,1),cell(nApp,1));
    [pF.xLabel.ind,pF.yLabel.ind] = deal(NaN);    
    [pF.xLabel.String,pF.yLabel.String] = deal(xStr,yStr);     
    
    % determines if the y-limits need to be recalculated
    if isempty(yLim) && setYLim
        [yLim,recalcYLim] = deal([1e10,-1e10],true);
    else
        recalcYLim = false;
    end
        
    % creates the sub-plots for all the axis
    for i = 1:nApp
        % creates the subplot and the graph of the metric
        hAx(i) = createSubPlotAxes(hP,[m,n],i);
        hold(hAx(i),'on');
        if m*n == 1; axis(hAx(i),'on'); end
        
        % creates the metric graph
        switch pType
            case ('polar')                
                % case is a polar plot
                
                % sets the patch colour index (first group only)
                if i == 1; iCol = mod(1:length(phi),2)+1; end
                
                % creates the polar plot axis
                polarPlotSetup(hAx(i),pF,false);
                
                % creates the polar plot patches
                indP = num2cell(1:length(phi));
                hPol = arrayfun(@(x,y,z)(createPolarPlotPatch...
                    (hAx(i),x,y,dphi,z)),phi,Y{i},iCol,'un',0);
                cellfun(@(h,i)(set(h,'UserData',i)),hPol,indP);
                
                % resets the x-axis limits
                xLim = get(hAx(i),'xlim');                                
            
            case ('boxplot')       
                % case is a boxplot
                if pP.plotErr
                    hPlt{i} = boxplot(Y{i},'sym','r*','outliersize',3);
                else
                    hPlt{i} = boxplot(Y{i},'sym','r');
                end
                
            case ('bar')
                % case is a bar graph
                hPlt{i} = bar(X,Y{i},'tag','hBar');
            
            case ('barstack')
                % case is a stacked bar graph
                hPlt{i} = bar(X,Y{i},'stacked');
        end
                
        % plots the fitted line (if required)        
        if pP.plotGrid; grid on; end               
        if ~isempty(Yfit);plot(X,Yfit{i},'r','linewidth',2); end
        set(hAx(i),'xlim',xLim,'UserData',i)
        
        % adds the error bars to the figure (if they exists & are required)
        if pP.plotErr && ~isempty(Ysem)
            if ~isempty(Ysem{i})
                addBarError(hPlt{i},X,Y{i},Ysem{i},'g'); 
            end
        end                   
        
        % sets the overall limits
        if recalcYLim
            yLimNw = get(hAx(i),'ylim');
            yLim = [min(yLim(1),yLimNw(1)),max(yLim(2),yLimNw(2))];                
        end        
        
        % formats the plot axis
        formatPlotAxis(hAx(i),pF,ind(i));                                   
    end
    
    % updates the axis limits
    if setYLim
        cellfun(@(x)(set(x,'ylim',yLim)),num2cell(hAx)); 
        if log10(diff(yLim)) > 3
            yLim = [max(1,yLim(1)),10^ceil(log10(yLim(2)))];
            cellfun(@(x)(set(x,'yscale','log','ylim',yLim)),num2cell(hAx));
        end
    end    
    
    % sets the non-aligned x/y labels
    formatMultiXYLabels(hAx,pF,[m,n]);
        
    % sets the legend (if first subplot and strings have been set)
    if ~isempty(lStr)
        % creates legend object
        [pF.Legend.String,pF.Legend.lgHorz] = deal(lStr,lgHorz);
        hLg = createLegendObj(hPlt{nApp},pF.Legend);        
        set(hLg,'Units','Normalized');
    end           
    
    % resets the axis positions
    resetAxesPos(hAx,m,n);           
    
    % sets the x-axis strings/ticklabels
    for i = 1:nApp
        % sets the axis properties        
        if rotXAx
            setGroupString(hAx(i),pF,X,xTickLbl,-90,0.01); 
        else
            set(hAx(i),'xTick',xTick,'xTickLabel',xTickLbl)
        end          
    end
        
    % resets the legend position
    if strcmp(pP.pMet,'Post-Contact (Movement Type)')                       
        % retrieves legend position and resets to the new position
        pLg = get(hLg,'position');
        set(hLg,'position',[(1-pLg(3))/2,1-pLg(4),pLg(3:4)])
                
        % updates the position of the legend object        
        resetLegendPos(hLg,hAx)
    elseif strcmp(pP.pMet,'Approach Angle (Polar Plot)')
        % creates the mean arrow for each sub-group
        for i = 1:nApp
            % converts the axis local coordinates to global coordinates
            axP = get(hAx(i),'position');
            
            % calculates the mean angle/vector length
            [phiMu,YMu] = deal(sum(Y{i}.*phi)/sum(Y{i}),mean(Y{i}));

            % creates the arrow annotation
            [x0,y0] = deal(axP(1)+axP(3)/2,axP(2));
            xArr = x0 + [0 axP(3)*YMu*cos((90-phiMu)*d2r)];
            yArr = y0 + [0 axP(4)*YMu*sin((90-phiMu)*d2r)];        

            % creates the arrow annotations
            hArr = annotation(hP,'arrow',xArr,yArr,...
                                 'color','r','linewidth',2);
            set(hArr,'units','pixels');
        end        
    end    
    
    % adds the fitted value strings (if the exists & are required)
    if ~isempty(fitVal) && pP.showPara
        % sets the font ratio
        sPPI = get(0,'ScreenPixelsPerInch');
        fR = max(0.8,min(newSz(3:4)./regSz(3:4))*sPPI/72);
        
        % sets the fitted strings
        for i = 1:nApp
            % sets focus to the corresponding subplot
            subplot(hAx(i))
            fSzT = setMinFontSize((14-2*szMx)*fR,'text');
            
            % retrieves the fit string text
            fitStr = getFittedParaText(fitVal{i});
            hText = text(0,0,fitStr,'FontWeight','bold','FontUnits',...
                                    'Pixels','FontSize',fSzT,'tag','Other');
            [hExt,dX,dY] = deal(get(hText,'Extent'),diff(xLim),diff(yLim));
            
            % calculates the new position of the parameters
            switch fitLoc
                case ('NorthEast')
                    hPosNw = [xLim(2)-hExt(3),yLim(2)-hExt(4)/2+pDel*dY,1];
%                     hPosNw = [(mean(xLim)-hExt(3)/2),(yLim(2)-hExt(4)),1];
                
                case ('SouthWest')
                    hPosNw = [xLim(1)+pDel*dX,yLim(1)+hExt(4)/2+pDel*dY,1];
            end
            
            % resets the textbox position
            set(hText,'position',hPosNw,'Units','Normalized')
        end
    end    
else
    % creates the subplot and the graph of the metric
    hAx = createSubPlotAxes(hP); 
    hold(hAx,'on'); 
    axis(hAx,'on');
        
    % creates the metric graph 
    switch pType
        case ('boxplot')       
            % case is a boxplot
            uData0 = get(hAx,'UserData');
            if pP.plotErr
                % case is plotting outliers
                hPlt = boxplot(Y,'sym','r*','outliersize',3);
            else
                % case is otherwise
                hPlt = boxplot(Y,'sym','r');
            end    
            
            % removes the extraneous text markers
            hText = findall(hAx,'type','text');
            if ~isempty(hText)
                if length(hText) == 1
                    tStr = {get(hText,'string')};
                else
                    tStr = get(hText,'String');
                end

                % sets the axis properties and deletes the text markers
                set(hAx,'UserData',uData0)
                delete(hText(~cellfun('isempty',tStr)));
            end
            
        case ('plot')
            % case is a line plot
            iPlt = num2cell(1:size(Y,2));
            hPlt = arr2vec(cellfun(@(y,i)...
                (plot(X,y,'UserData',i)),num2cell(Y,1),iPlt));
            
        case ('bar')
            % case is a bar graph
            hPlt = bar(X,Y,'tag','hBar');
            
        case ('barstack')
            % case is a stacked bar graph
            hPlt = bar(X,Y,'stacked');            
    end
       
    % adds in the plot errorbars (if they exist and are required)
    if pP.plotErr && ~isempty(Ysem)
        addBarError(hPlt,X,Y,Ysem,'g'); 
    end                     
    
    % sets the axis properties
    set(hAx,'Units','Normalized','xLim',xLim,'UserData',1);
    if pP.plotGrid; grid on; end 
    if ~isempty(yLim); set(hAx,'ylim',yLim); end
    
    yLim = get(hAx,'yLim');
    if log10(diff(yLim)) > 3
        yLim = [max(1,yLim(1)),10^ceil(log10(yLim(2)))];
        set(hAx,'yscale','log','ylim',yLim);
    end     
    
    % formats the plot axis
    pF.Title(1).Font.FontSize = 30;
    [pF.xLabel.String,pF.yLabel.String] = deal(xStr,yStr);
    formatPlotAxis(hAx,pF,1);     
    
    % resets the axis positions
    resetAxesPos(hAx,1,1);    
    
    % special additions to the figure
    switch pP.pMet
        case ('Approach Angle (CDF)')
            % case is the Approach Angle (CDF)
            plot([0 0],get(hAx,'ylim'),'r--','HitTest','off')
            plot(get(hAx,'xlim'),50*[1 1],'r--','HitTest','off')

        case ('Approach Angle (Vector)')
            % case is the Approach Angle (Vector)
            set(hAx,'ylim',[-0.05 1.05])
    end     
    
    % sets up the legend (if the subplot has one)
    if ~isempty(lStr)
        % creates the legend object
        [pF.Legend.String,pF.Legend.lgHorz] = deal(lStr,lgHorz);
        hLg = createLegendObj(hPlt,pF.Legend);  
        set(hLg,'Units','Normalized');        
        
        % specialised legend placement
        pLg = get(hLg,'position');
        switch pP.pMet
            case ('Approach Angle (CDF)')     
                % resets the legend position            
                set(hLg,'Position',[(1-pLg(3)),0.5*(1-pLg(4)),pLg(3:4)])   

                % updates the position of the legend object        
                resetLegendPos(hLg,hAx)
            case ('Post-Contact (Movement Type)')                
                % retrieves legend position and resets to the new position
                set(hLg,'position',[(1-pLg(3))/2,1-pLg(4),pLg(3:4)])
                
                % updates the position of the legend object        
                resetLegendPos(hLg,hAx)                
                
            otherwise
                % sets the fixed location of the legend
                set(hLg,'Location',lgLoc)    
        end            
    end      
    
    % sets the x-axis text labels
    if rotXAx
        setGroupString(hAx,pF,X,xTickLbl,30,0.0); 
    elseif setXTick
        set(hAx,'xTick',xTick,'xTickLabel',xTickLbl)
    end                             
    
    % resets the legend position
    if strcmp(pP.pMet,'Post-Contact (Movement Type)')
        % retrieves legend position and resets to the new position
        pLg = get(hLg,'position');
        [pAx,dH] = deal(get(hAx,'position'),(lgOfs+pLg(4)));
        
        % updates the axis dimensions for all subplots
        pAx(4) = pAx(4) - dH;
        set(hAx,'position',pAx)

        % updates the legend position
        pLgNw = [(0.5-pLg(3)/2),1-(lgOfs+pLg(4)),pLg(3:4)];                
        set(hLg,'position',pLgNw)
    end                                               
end
    
% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% resets the mean/sem shape/slope factors
for i = 1:length(plotD)
    [plotD(i).ATD_mn,plotD(i).ATD_sem] = ...
                                deal({plotD(i).ATD(1)},{plotD(i).ATD(2)});
    [plotD(i).nTD_mn,plotD(i).nTD_sem] = ...
                                deal({plotD(i).nTD(1)},{plotD(i).nTD(2)});
                            
    plotD(i).TM = cellfun(@(x)(x(:)),plotD(i).TM,'un',0);
    plotD(i).TE = cellfun(@(x)(x(:)),plotD(i).TE,'un',0);
    plotD(i).DE = cellfun(@(x)(x(:)),plotD(i).DE,'un',0);
end

% ----------------------------------------------------------------------- %
% ---                         OTHER FUNCTIONS                         --- %
% ----------------------------------------------------------------------- %

% --- calculates the approach angle of the flies (wrt to the circle centre)
function PhiA = calcApproachAngle(xB,yB)

% calculates the angle from final point in the sequence to circle centre
PhiF = atan2(yB(end),xB(end));

% calculates the rotated coordinates
[xBR,yBR] = rotateCoords(xB,yB,PhiF);
xBR = xBR - xBR(end);

% calculates the mid-distance of each line
D = sqrt((0.5*(xBR(2:end)+xBR(1:end-1))).^2+(0.5*(yBR(2:end)+yBR(1:end-1))).^2);

% determines the last frame in the sequence where the difference in rotated
% x-coordinates is less than zero. if there are none, then set the index to
% the first frame (i.e., all frames are included in the calculation)
i0 = find(diff(xBR) < 0,1,'last');
if isempty(i0); i0 = 0; end

% otherwise, calculate the mean approcach angle from the sequence
if i0 == length(xB)
    PhiA = NaN;
else
    ii = (i0+1):length(xB);
    Q = 1./D(ii(1:(end-1))); Q = Q/sum(Q);
    [dyBR,dxBR] = deal(diff(yBR(ii)),diff(xBR(ii)));    
    PhiA = sum(Q.*atan2(dyBR,dxBR).*(abs(dyBR) > 1e-6));
end

% --- calculates the Approach Angle (CDF) values
function [xCDF,fCDF] = calcAngleCDF(PhiA)

% initialisations
xCDF = (-90:0.1:90)';

% calculates the empirical CDF 
if isempty(PhiA)
    fCDF = NaN(size(xCDF));
else
    % determines the unique values from the x array and reset the values
    [f,x] = ecdf(PhiA);
    [x,iA] = unique(x); f = f(iA);

    % determines the correct values
    ii = find(x >= -90,1,'first'):find(x <= 90,1,'last');

    % reshapes the array so that the values 
    [f,x] = deal((f(ii) - f(ii(1)))/(f(ii(end)) - f(ii(1))),x(ii));

    % set the interpolation array
    if length(f) < 2
        fCDF = double(xCDF > x);
    else
        fCDF = min(1,max(0,interp1(x,f,xCDF,'linear','extrap')));
    end
end

% --- sets up the fitted parameter text string
function fitStr = getFittedParaText(p)

% initialisations
latStr = {'mu','sigma'};
[fitStr,fStr] = deal([],fieldnames(p));

% sets the strings for all the parameters
for i = 1:length(fStr)
    % retrieves the new value
    nwVal = eval(sprintf('p.%s',fStr{i}));
    
    % sets the values into the overall string
    if length(nwVal) == 1
        fitStr = sprintf('%sR^2 = %.4f',fitStr,nwVal(1));
    else
        if any(strcmp(latStr,fStr{i}))
            if isnan(nwVal(2))
                fitStr = sprintf('%s\\%s = %.2f',fitStr,fStr{i},nwVal(1));                
            else
                fitStr = sprintf('%s\\%s = %.2f +/- %.2f',...
                                        fitStr,fStr{i},nwVal(1),nwVal(2));
            end
        else
            if isnan(nwVal(2))
                fitStr = sprintf('%s%s = %.2f',fitStr,fStr{i},nwVal(1));                
            else            
                fitStr = sprintf('%s%s = %.2f +/- %.2f',...
                                        fitStr,fStr{i},nwVal(1),nwVal(2));
            end
        end
    end
    
    % adds a carriage return (except for last parameter)
    if i ~= length(fStr); fitStr = sprintf('%s\n',fitStr); end
end
