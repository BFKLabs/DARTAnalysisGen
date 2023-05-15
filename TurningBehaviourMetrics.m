% --- calculates the fly turning behaviour metrics over the duration
%     of an experiment (2D experiment only)
function pData = TurningBehaviourMetrics(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Turning Behaviour (Metrics)';
pData.Type = {'Pop','Multi'};
pData.fType = [2 1 1 3];
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
rI.Shape = '2D (Circle)';
rI.Stim = 'None';
rI.Spec = 'None';

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 4;                       
cP = setParaFields(nPara);

% sets the parameter fields
cP(1) = setParaFields([],'Number',10,'nFrm','Turn Duration (Frames)',[7 15 true]);
cP(2) = setParaFields([],'Number',0.2,'vTol','Speed Threshold (mm/s)',[0.01 2.00 false]);
cP(3) = setParaFields([],'Number',5,'dTol','Distance Threshold (mm)',[1 25 false]);
cP(4) = setParaFields([],'Number',2,'dR','Outside Edge Distance (mm)',[0.1 5 false]);

% sets the tool-tip strings
cP(1).TTstr = 'Number of frames over which a turning event is determined';
cP(2).TTstr = 'Minimum speed the fly must have after coming in contact with wall';
cP(3).TTstr = 'Minimum distance to be travelled over a turning event';
cP(4).TTstr = 'The distance from outside edge whereby fly is considered in contact with wall';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 1;
pP = setParaFields(nPara);

% sets the parameter fields
pP(1) = setParaFields([],'Boolean',0,'plotGrid','Show Axis Gridlines');

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
nSub = 2;
pF = setFormatFields(nSub);

% initialises the font structs
pF.Title = setFormatFields(setupFontStruct('FontSize',24),'',nSub);
pF.xLabel = setFormatFields(setupFontStruct('FontSize',20),' ',nSub);
pF.yLabel = setFormatFields(setupFontStruct('FontSize',20),'',nSub);
pF.Axis = setFormatFields(setupFontStruct('FontSize',16),[]);

% sets the total turn subplot titles
pF.Title(1).String = 'Total Turns';
pF.yLabel(1).String = 'Turn Count';

% sets the direction changes subplot titles
pF.Title(2).String = 'Direction Changes';
pF.yLabel(2).String = '% Changes';

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);
[Type1,Type2,Stats] = deal(3,2,{'Comp'});

% sets the dependent output variables
oP = addYVarField(oP,'Turn Count','NT',Stats,Type1,[],1);
oP = addYVarField(oP,'Turn Count (Mean)','NT_mn',[],Type2);
oP = addYVarField(oP,'Turn Count (SEM)','NT_sem',[],Type2);
oP = addYVarField(oP,'Turn Count (Median)','NT_md',[],Type2);
oP = addYVarField(oP,'Turn Count (Lower Quartile)','NT_lq',[],Type2);
oP = addYVarField(oP,'Turn Count (Upper Quartile)','NT_uq',[],Type2);
oP = addYVarField(oP,'Direction Change','dDC',Stats,Type1,[],1);
oP = addYVarField(oP,'Direction Change (Mean)','dDC_mn',[],Type2);
oP = addYVarField(oP,'Direction Change (SEM)','dDC_sem',[],Type2);
oP = addYVarField(oP,'Direction Change (Median)','dDC_md',[],Type2);
oP = addYVarField(oP,'Direction Change (Lower Quartile)','dDC_lq',[],Type2);
oP = addYVarField(oP,'Direction Change (Upper Quartile)','dDC_uq',[],Type2);
             
% --- sets the data cursor update function
function dTxt = dataCursorFunc(hObj,evnt,dcObj)

% updates the plot data fields
dcObj.getCurrentPlotData(evnt);

% retrieves the current plot data
sP = retParaStruct(dcObj.pData.sP);

% field retrievals
uStr = {'count','%'};
mStr = {'Turn Count','Direction % Changes'};
iAx = dcObj.getSelectAxesIndex;
grpName = dcObj.pData.appName(sP.Sub.isPlot);

% sets the common class fields
dcObj.pType = 'Boxplot';
dcObj.yName = mStr{iAx};
dcObj.yUnits = uStr{iAx};
dcObj.grpName = grpName;
[dcObj.useGrpHdr,dcObj.combFig] = deal(false);
[dcObj.xName,dcObj.xGrp] = deal('Group Name',grpName);

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

% initialises the plot value data struct
plotD = initPlotValueStruct(snTot,pData,cP,...
                                 'NT',[],'NT_mn',[],'NT_sem',[],...
                                 'NT_md',[],'NT_lq',[],'NT_uq',[],...
                                 'dDC',[],'dDC_mn',[],'dDC_sem',[],...
                                 'dDC_md',[],'dDC_lq',[],'dDC_uq',[]);                        

% other memory allocations                       
tInd = cell(nExp,nApp);

% parameters
nS1 = 2:4;
[nSm,nS2] = deal(5,(nS1(end)+1):cP.nFrm);

% creates the waitbar figure
wStr = {'Turn Detection','Turning Metric Calculations'};
h = ProgBar(wStr,'Turning Behaviour Calculations');

% ------------------------------------------------------- %
% --- INTER-FRAME DISTANCE CALCULATION & THRESHOLDING --- %
% ------------------------------------------------------- %
    
% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure 
    wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
    if h.Update(1,wStrNw,i/(1+nExp))
        [plotD,ok] = deal([],false);
        return
    end

    % calculates the video frame rate and experiment apparatus indices  
    sFac = snTot(i).sgP.sFac;
    iApp = find(~cellfun('isempty',snTot(i).iMov.flyok));
    nFrm = size(snTot(i).Px{iApp(1)},1);
    Rtol = cP.vTol/sFac;
    
    % sets the relevant x/y-locations for the current experiment  
    [dPx0,dPy0,R0] = get2DCoordsBG(snTot(i),iApp);
    if ~iscell(dPx0)
        [dPx0,dPy0,R0] = deal({dPx0},{dPy0},{R0});
    end
    
    % calculates the scaled coordinates/radii
    dPx = cellfun(@(x)(x/sFac),dPx0,'un',0);
    dPy = cellfun(@(x)(x/sFac),dPy0,'un',0);
    R = cellfun(@(x)(x/sFac),R0,'un',0);
    
    % determines the turning events for each frame
    for j = 1:length(iApp)
        % calculates the angle of the fly wrt the circle centre
        tInd{i,iApp(j)} = cell(1,size(dPx{j},2));        
        for k = 1:size(dPx{j},2)
            kGrp = getGroupIndex(~isnan(dPx{j}(:,k)));
            for kk = 1:length(kGrp)
                dPx{j}(kGrp{kk},k) = smooth(dPx{j}(kGrp{kk},k),nSm);
                dPy{j}(kGrp{kk},k) = smooth(dPy{j}(kGrp{kk},k),nSm);
            end
        end
        
        % calculates the inter-frame displacement/angle 
        [d2Px,d2Py] = deal(diff(dPx{j},[],1),diff(dPy{j},[],1));
        dD = [zeros(1,size(dPx{j},2));sqrt(d2Px.^2 + d2Py.^2)];
        
        % calculates the distance from the circle centre, and from this
        % determines which flies are considered in contact with the edge
        Rho = sqrt(dPx{j}.^2 + dPy{j}.^2);
        onEdge = Rho > (repmat(R{j}(:)',nFrm,1) - cP.dR/sFac);
                
        % determines the index groups where the fly is in contact with the
        % edge of the arena. groupings 
        iGrp = cellfun(@(x)(getGroupIndex(x)),num2cell(onEdge,1),'un',0);                
        dDCnw = NaN(1,length(iGrp));
        for k = 1:length(iGrp)
            % memory allocations
            if ~isempty(iGrp{k})
                % determines the first index of each group
                ind1 = cellfun(@(x)(x(1)),iGrp{k});
                ind1 = ind1(ind1 > 1);
                
                % determines which groups are when there is a significant
                % movement towards the outside (must be greater than vTol)
                if ~isempty(ind1)
                    ii = cellfun(@(x)...
                            (diff(Rho(x+(-1:0),k))>Rtol),num2cell(ind1));
                    [ind1,iGrp{k}] = deal(ind1(ii),iGrp{k}(ii));
                end
                
                % removes the groups that either within the frame tolerance
                % of the previous group
                ii = ind1 < (nFrm - cP.nFrm);
                
                % removes any groups that have start indices within the frame
                % tolerance of the preceding group
                if any(ii)
                    iPr = find(ii,1,'first');                                
                    for iNw = (iPr+1):length(iGrp{k})
                        % only check valid groups
                        if ii(iNw)
                            % determines if the initial index of the new group
                            % relative to the previous is greater than the
                            % frame tolerance
                            if ((ind1(iNw) - ind1(iPr)) > cP.nFrm)
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
                    [ind0,iGrp{k}] = deal(num2cell(ind1(ii)-1),iGrp{k}(ii));

                    % from the remaining groups, determine which events have the
                    % fly travel greater than the distance tolerance (dTol)
                    if ~isempty(ind0)
                        % calculates the distance over the turn path and
                        % ensures that distance meets the tolerance
                        DGrp = cellfun(@(x)(sum(dD(x+(0:cP.nFrm),k))),ind0);                                                
                        jj = DGrp > cP.dTol/sFac;
                        
                        % reorders the 
                        tInd{i,iApp(j)}{k} = cell2mat(ind0(jj));          
                        iGrp{k} = iGrp{k}(jj);
                    end
                end            
            end
                        
            % if any turns have been found, then calculate the associated
            % metrics with the turns 
            if ~isempty(tInd{i,iApp(j)}{k})
                % sets the start indices of the turn events, and the
                % indices of the subsequent turn paths (nFrm frames long)
                iNw = num2cell(tInd{i,iApp(j)}{k});               
                jNw = cellfun(@(x)(x+(0:cP.nFrm)),iNw,'un',0);
                
                % calculates the flies orientation (maps the angles to the
                % domain [0 2*pi]
                dPxP = cell2mat(cellfun(@(x)(dPx{j}(x,k)),jNw,'un',0)');
                dPyP = cell2mat(cellfun(@(x)(dPy{j}(x,k)),jNw,'un',0)');
                Phi = calcFlyOrientation(dPxP,dPyP);
                
                % calculates the entry direction (positive for left heading
                % entry, negative for right heading entry)
                Phi0 = atan2(cellfun(@(x)(dPy{j}(x+1,k)),iNw),...
                             cellfun(@(x)(dPx{j}(x+1,k)),iNw))';
                dEntry = sign(atan2(sin(Phi(1,:)-Phi0),cos(Phi(1,:)-Phi0)));
            
                % calculate the turning paths angle wrt the initial angle
                Q = Phi - repmat(Phi(1,:),size(Phi,1),1);
                dPhi = unwrap(atan2(sin(Q),cos(Q)),[],1);
                dPhi = repmat(dEntry,size(dPhi,1),1).*dPhi;                

                % ensures the turn angles are within the [-pi pi] range
                dQ = dPhi(2:end,:) - dPhi(1:end-1,:);
                d2Phi = abs(atan2(sin(dQ),cos(dQ)));                
                for jj = find(any(d2Phi > pi/2,1))
                    for i0 = find(d2Phi(:,jj)>pi/2,1,'first'):size(d2Phi,1)
                        % determines if there is a large jump between
                        % frames
                        if d2Phi(i0,jj) > pi/2
                            % if so, remove the jump
                            if (dPhi(i0+1,jj) > 0)
                                dPhi(i0+1,jj) = dPhi(i0+1,jj) - pi;
                            else
                                dPhi(i0+1,jj) = dPhi(i0+1,jj) + pi;
                            end
                            
                            % updates the next 
                            if (i0 < size(d2Phi,1))
                                dQ = dPhi(i0+2,jj) - dPhi(i0+1,jj);
                                d2Phi(i0+1,jj) = abs(atan2(sin(dQ),cos(dQ)));
                            end
                        end
                    end
                end
                                                                
                % determines which paths have changed their direction. from
                % this calculate the proportion change in direction
                sdPhi = sign(unwrap(dPhi,1));
                isDC = (sum(sdPhi(nS1,:),1).*sum(sdPhi(nS2,:),1)) < 0;
                dDCnw(k) = sum(isDC)/size(Phi,2);
            end
        end
        
        % appends the overall turn count to the total array
        ifok = snTot(i).iMov.flyok{iApp(j)};
        plotD(iApp(j)).NT(1,ifok,i) = num2cell(cellfun('length',iGrp)); 
        plotD(iApp(j)).dDC(1,ifok,i) = num2cell(dDCnw);        
    end    
end
    
% sets/calculates the raw, mean and SEM proportional activity values
for i = 1:nApp
    % updates the waitbar figure (if more than one solution file)
    wStrNw = sprintf('%s (Apparatus %i of %i)',wStr{2},i,nApp);
    if h.Update(2,wStrNw,i/(1+nApp))
        [plotD,ok] = deal([],false);
        return
    end    
    
    % calculates the metric statistics
    plotD(i) = calcMetricStats(plotD(i),'NT'); 
    plotD(i) = calcMetricStats(plotD(i),'dDC');
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

% memory allocation
nSub = 2;
[hAx,xi] = deal(cell(1,nSub),1:length(ind));

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the formatting struct based on the plot type
axSz = 16 - max(0,nApp-4);  

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %
   
% creates the subplots for each metric type
for i = 1:2
    % sets the plot values
    switch (i)
        case (1) % case is the number of turns
            Yplt = combineNumericCells(field2cell(p,'NT'));
        case (2) % case is the direction change proportion
            Yplt = 100*combineNumericCells(field2cell(p,'dDC'));
    end    
    
    % sets focus to the appropriate axes
    hAx{i} = createSubPlotAxes(hP,[1,nSub],i);
    
    % creates the boxplot and removes all the text labels
    boxplot(Yplt);
    delete(findall(hAx{i},'type','text'))
    
    % sets the axis properties
    xLim = 0.5 + [0,size(Yplt,2)];
    set(hAx{1},'xtick',xi,'TickLength',[0 0])
    if (i == 2)
        set(hAx{i},'ylim',[-2 103],'xLim',xLim)
    else
        set(hAx{i},'ylim',[0 max(get(hAx{i},'ylim'))],'xLim',xLim)
    end
    
    % sets the x-axis labels
    set(hAx{i},'UserData',i);
    formatPlotAxis(hAx{i},pF,i);              
    
    % adds in the gridlines (if checked)    
    if pP.plotGrid; grid(hAx{i},'on'); end    
end

% formats and resets the axis positions
resetAxesPos(hAx,1,2); 

% sets the group strings
for i = 1:2
    setGroupString(hAx{i},pF,xi,snTot(1).iMov.pInfo.gName(ind),-90);
end
    
% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)
