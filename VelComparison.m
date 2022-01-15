% --- comparison of the average pre/post-stimuli speeds for different time 
%     throughout the day (long experiment only)
function pData = VelComparison(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Pre & Post Stimuli Speed Comparison';
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
    [pData.hasSR,pData.hasRS] = deal(true,false);
    pData.sP = initSpecialPara(snTotL,pData,'nGrp',0);  
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
nPara = 3;                       
cP = setParaFields(nPara);
isLong = ~detIfShortExpt(field2cell(snTot,'T'));

% sets the tab list names
a = {'1 - General'};

% sets the parameter list for the drop-down list
pList = cellfun(@num2str,num2cell([1 2 4 6 8 12]),'un',0);

% sets the parameter fields
cP(1) = setParaFields(a{1},'List',{1+2*isLong,pList},'nGrp','Number of Daily Time Groups',[],{0,~isLong});
cP(2) = setParaFields(a{1},'Number',60,'tBefore','Pre-Stimuli Duration (sec)',[10 300 false]);
cP(3) = setParaFields(a{1},'Number',60,'tAfter','Post-Stimuli Duration (sec)',[10 300 false]);
cP(4) = setParaFields(a{1},'Boolean',0,'useZG','Use Zeitgeiber For Group Times');

% sets the other tool-tip strings
cP(1).TTstr = 'The number of groups that the day is split up into';
cP(2).TTstr = 'Time period before stimuli used to calculate average speed';
cP(3).TTstr = 'Time period after stimuli used to calculate average speed';

% adds the unique motor parameters
cP = addUniqueMotorPara(cP,snTot);

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = 9;                       
pP = setParaFields(nPara);

% sets the parameter list strings
pList = {'Bar Graph','Boxplot'};

% sets the tab list names
a = {'1 - Scatterplot','2 - Metrics'};

% sets the parameter fields
pP(1) = setParaFields(a{2},'List',{1,pList},'pType','Plot Type');
pP(2) = setParaFields(a{1},'Number',6,'mSize','Plot Marker Size',[1 20 1]);
pP(3) = setParaFields(a{2},'Number',0.8,'pW','Bar Graph Relative Width',[0 1 false],{1,1});
pP(4) = setParaFields(a{2},'Boolean',1,'grpTime','Group Data By Time Group');
pP(5) = setParaFields(a{2},'Boolean',1,'plotErr','Show Error Bars/Outliers');
pP(6) = setParaFields(a{1},'Boolean',0,'plotPolar','Plot Polar Coordinates');
pP(7) = setParaFields(a{1},'Boolean',1,'addTrend','Add Regression Line & Equation',[],{6,1});
pP(8) = setParaFields(a{2},'Boolean',0,'plotGrid','Plot Bar Graph/Boxplot Gridlines');
pP(9) = setParaFields(a{2},'Boolean',1,'plotFixedY','Fixed Subplot Limits To Overall Maximum');

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
[nApp,nSub] = deal(length(snTot.iMov.ok),3);
pF = setFormatFields(nSub);

% initialises the font structs
pF.Title = setFormatFields(setupFontStruct('FontSize',20),'',nSub);
pF.xLabel = setFormatFields(setupFontStruct('FontSize',14),'',nSub);
pF.yLabel = setFormatFields(setupFontStruct('FontSize',14),'',nSub);
pF.Axis = setFormatFields(setupFontStruct('FontSize',12),[]);

% sets x/y labels strings for the 1st plot
pF.xLabel(1).String = 'Pre-Stimuli Speed (mm s^{-1})';
pF.yLabel(1).String = 'Post-Stimuli Speed (mm s^{-1})';

% sets x/y labels strings for the 2nd plot
pF.yLabel(2).String = 'Speed (mm s^{-1})';

% sets x/y labels strings for the 2nd plot
pF.yLabel(3).String = 'Speed (mm s^{-1})';

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Axis(1).String{i} = snTot.iMov.pInfo.gName{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot);
[Stats,Type] = deal({'CompMulti','Tgrp'},[1 3]);

% sets the independent output variables
oP = addXVarField(oP,'Time Group','Tgrp','Group');

% sets the dependent output variables
oP = addYVarField(oP,'Pre-Stim Spd','V1',Stats,Type,{'Tgrp'},1);
oP = addYVarField(oP,'Post-Stim Spd','V2',Stats,Type,{'Tgrp'},1);
oP = addYVarField(oP,'Speed Ratio','M',Stats,Type,{'Tgrp'},1);

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
cP.movType = 'Absolute Distance';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensioning and memory allocation
nGrp = str2double(cP.nGrp);
[nApp,nExp,ok] = deal(length(snTot(1).iMov.ok),length(snTot),true);

% fixed parameters
[tBefore,tAfter] = deal(cP.tBefore,cP.tAfter);

% sets the daily time group strings
Tgrp = setTimeGroupStrings(nGrp,cP.Tgrp0,cP.useZG); 

% retrieves the other calculation parameters (if they exist)
[devType,chType] = deal([]);
if isfield(cP,'devType'); devType = cP.devType; end
if isfield(cP,'chType'); chType = cP.chType; end

% memory allocation
plotD = initPlotValueStruct(snTot,pData,cP,NaN(1,nGrp),'Tgrp',Tgrp,...
                                'V1',[],'V1_mn',[],'V1_sem',[],...
                                'V2',[],'V2_mn',[],'V2_sem',[],...
                                'M',[],'M_mn',[],'M_sem',[]);                                     
                             
% ---------------------------------------------- %
% --- PRE/POST STIMULI VELOCITY CALCULATIONS --- %
% ---------------------------------------------- %

% sets the waitbar offset (is > 0 for more than one 
[wStr,wOfs] = deal({'Overall Progress'},1);
if nargin == 5
    % resets the strings
    [h,wOfs] = deal(varargin{1},wOfs+1);    
    wStr = [h.wStr(1),wStr];   
    
    % resets the waitbar figure fields
    for i = 2:length(wStr); h.Update(i,wStr{i},0); end

    % collapses the waitbar figure
    hh0 = getHandleSnapshot(findall(h));
    h.collapseProgBar(1);   
    
else
    % creates the waitbar figure
    h = ProgBar(wStr,'Pre/Post Stimuli Velocity Calculations');    
end

% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure (if more than one solution file)
    wStrNw = sprintf('%s (Experiment %i of %i)',wStr{wOfs},i,nExp);
    if h.Update(wOfs,wStrNw,i/(1+nExp))
        [plotD,ok] = deal([],false);
        return
    end
        
    % sets the video/stimuli time stamps into a single vector
    flyok = snTot(i).iMov.flyok;
    Ttot = cell2mat(snTot(i).T);
    Ts = getMotorFiringTimes(snTot(i).stimP,devType,chType);
    if ~isempty(Ts)
        % determines the indices of the stimuli events within the total
        % experiment, and determines what time groups that the stimuli
        % events took place in         
        T0 = snTot(i).iExpt(1).Timing.T0;
        indGrp = detTimeGroupIndices(Ts,T0,nGrp,cP.Tgrp0,true);
        
        % calculates the pre-stimuli velocities
        indV1 = cellfun(@(x)(find(Ttot>(x-tBefore),1,'first'):...
                    find(Ttot<x,1,'last')),num2cell(Ts),'un',0);
        indV2 = cellfun(@(x)(find(Ttot>x,1,'first'):find(...
                    Ttot<(x+tAfter),1,'last')),num2cell(Ts),'un',0);                    
                
        % calculates the pre/post stimuli velocities for all flies, and
        % bins the values according to their time group bin
        for j = 1:nApp                
            % calculates the pre-stimuli velocities   
            if (~isempty(snTot(i).Px{j}))
                V1 = calcBinnedFlyMovement(snTot(i),Ttot,indV1,cP,j,flyok{j},1);
                V1 = cell2mat(cellfun(@(x)(x/tBefore),V1,'un',0));
                                
                % calculates the post-stimuli velocities                    
                V2 = calcBinnedFlyMovement(snTot(i),Ttot,indV2,cP,j,flyok{j},1);
                V2 = cell2mat(cellfun(@(x)(x/tAfter),V2,'un',0));
                                
                % sets the raw data values
                plotD(j) = setRawDataValues(plotD(j),snTot(i),V1,indGrp,'V1',i,j,3);
                plotD(j) = setRawDataValues(plotD(j),snTot(i),V2,indGrp,'V2',i,j,3);                                                      
            end
        end
    end    
end

% sets the final pre/post stimuli event velocities (for each apparatus)
for i = 1:nApp            
    if (~all(cellfun(@isempty,plotD(i).V1(:))))
        % sets/calculates the plot values
        plotD(i).M = cellfun(@(x,y)((x./y)),plotD(i).V2,plotD(i).V1,'un',0);
        
        % calculates the metric statistics for each of the values
        plotD(i) = calcMetricStats(plotD(i),'V1',3,false);        
        plotD(i) = calcMetricStats(plotD(i),'V2',3,false);
        plotD(i) = calcMetricStats(plotD(i),'M',3,false);                                  
    end    
end
    
% updates the waitbar figure
if ~h.Update(wOfs,'Pre/Post Stimuli Velocity Calculations Complete!',1)
    if nargin == 5
        resetHandleSnapshot(hh0)
    else
        h.closeProgBar()
    end
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

% if the current data set is empty, then exit the loop
if (all(isnan(plotD{1}(sP.pInd).V1_mn))); return; end

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the indices to be plotted
iPlotT = find(sP.pT);
if (isempty(iPlotT)); return; end

% sets the parameters from the data structs
p = plotD{1}(sP.pInd);

% sets the time group value and the x-plot values
nGrp = str2double(cP.nGrp);
xi = 1:nGrp;

% other parameters
isBar = strcmp(pP.pType,'Bar Graph');
[c,mm,lWid,nTick,pWL] = deal('rk','+o*xsd^v><ph',2,6,0.85);

% retrieves the panel object handle
hP = getCurrentAxesProp('Parent');

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% sets the legend strings
[lStr1,lStr2] = deal(p.Tgrp,{'Pre-Stimuli','Post-Stimuli'});   
pF.Legend(1).String = lStr1(iPlotT);

% sets the main title string
pF.Title(1).String = pF.Axis(1).String{sP.pInd};
a = (~cellfun(@isempty,p.V1) & ~cellfun(@isempty,p.V2));

% sets the polar plot labels
if (pP.plotPolar)
    pP.addTrend = false;
    pF.yLabel(1).String = 'Distance (mm s^{-1})';
    pF.xLabel(1).String = 'Phase Angle (deg)';
end

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% --- AVERAGE SPEED VS GROUP --- %
% ------------------------------ %

% sets up the sub-plot
hAx = createSubPlotAxes(hP,[2,2],4); 
hold(hAx,'on');
set(hAx,'Units','Normalized','box','on')

% creates the bar graph/boxplot figures
[hPlot,xTick] = plotDoubleBarBoxMetrics(hAx,p,{'V1','V2'},pP);      

% creates the subplot legend
if (pP.grpTime)
    % resets the legend/x-axis tick labels
    [pF.Legend.String,pF.Legend.lgHorz] = deal(lStr2,true);      
    if (isBar)
        % case is the graph is a bar graph
        [xLim,xi] = deal([1 nGrp],1:nGrp);
    else
        % case is the graph is a boxplot graph
        [xLim,xi] = deal(xTick([1 end]) + 0.5*[-1 1],xTick);
        set(hAx,'xTick',xTick)
    end    
    
    % updates the graph properties
    set(hAx,'xlim',xLim + 0.5*[-1 1])    
else
    % resets the legend/x-axis tick labels
    pF.Legend.String = lStr1;  
    if (isBar); xTick = [1 2]; end           
    
    % updates the graph properties
    set(hAx,'xTick',xTick,'xticklabel',lStr2)               
end

% updates the axis properties
formatPlotAxis(hAx,pF,2); 

% creates the legend object
hLg = createLegendObj(hPlot,pF.Legend);
    
% turns the grid on (if specified)
if (pP.plotGrid); set(hAx,'ygrid','on'); end

% sets the 
if (pP.plotFixedY)
    % determines the maximum y-axis extent over all apparatus
    if (isBar)
        V1mx = max(cellfun(@(x,y)(detOverallLimit(x+y)),field2cell(...
                            plotD{1},'V1_mn'),field2cell(plotD{1},'V1_sem')));
        V2mx = max(cellfun(@(x,y)(detOverallLimit(x+y)),field2cell(...
                            plotD{1},'V2_mn'),field2cell(plotD{1},'V2_sem')));
    else
        V1mx = max(cellfun(@(x)(detOverallLimit(x)),field2cell(plotD{1},'V1')));
        V2mx = max(cellfun(@(x)(detOverallLimit(x)),field2cell(plotD{1},'V2')));
    end
                        
    % sets the axis limits for the two speed comparison subplots
    setStandardYAxis(hAx,[],nTick,detOverallLimit([V1mx V2mx]));
end

% resets the axes position
resetAxesPos(hAx,1,1);

% resets the location of the legend object
[axP,lgP] = deal(get(hAx,'position'),get(hLg,'position'));
if (pP.grpTime)
    % case is that time is being grouped
    lgP(1) = axP(1) + 0.5*(axP(3)-lgP(3));
    lgP(2) = sum(axP([2 4]));
    set(hLg,'position',lgP);
else
    % case is that time is not being grouped
    lgP(2) = axP(2) + 0.5*(axP(4)-lgP(4));
    set(hLg,'position',lgP)    
    
    %
    if isBar
        hErr = findall(hAx,'tag','hErr');                
        cSzNw = hErr(1).CapSize*(axP(3)-pWL*lgP(3))/axP(3);
        
        for i = 1:length(hErr)
            hErr(i).CapSize = cSzNw;
        end
    end    
    
    % resets the position of the axis
    resetObjPos(hAx,'width',-lgP(3),1)    
end

% adds the x-axis labels 
if (pP.grpTime); setGroupString(hAx,pF,xi,lStr1,-90); end

% --- AVERAGE SPEED VS PRE/POST STIMULI --- %
% ----------------------------------------- %

% sets up the sub-plot
hAx2 = createSubPlotAxes(hP,[2,2],2); 
hold(hAx2,'on');

% creates the bar graph/boxplot
yLim = plotBarBoxMetrics(hAx2,1:nGrp,p,'M',pP,[]);             
    
% updates the axis properties
pF.yLabel(3).String = 'Gradient (Unitless)';
set(hAx2,'xtick',get(hAx,'xtick'),'xticklabel',[],...
         'xLim',[1 nGrp]+0.5*[-1 1],'box','on')        
    
% formats the axis plot
formatPlotAxis(hAx2,pF,3); 

% turns the grid on (if specified)
if (pP.plotGrid); set(hAx2,'ygrid','on'); end

% sets the 
if (pP.plotFixedY)
    % determines the maximum y-axis extent over all apparatus
    if (isBar)
        Mmx = max(cellfun(@(x,y)(detOverallLimit(x+y)),field2cell(...
                            plotD{1},'M_mn'),field2cell(plotD{1},'M_sem')));
    else
        Mmx = max(cellfun(@(x)(detOverallLimit(x)),field2cell(plotD{1},'M')));
    end                      
    
    % sets the axis limits for the two speed comparison subplots
    if (Mmx > 500)        
        set(hAx2,'ylim',[0 Mmx],'yscale','log');
    else    
        setStandardYAxis(hAx2,[],nTick,Mmx);        
    end
else
    % determines if the y-axis scale is extremely large
    if (max(get(hAx2,'ylim')) > 500)
        % if so, then use a log scale
        set(hAx2,'yscale','log');
    else
        % otherwise, use a standard axis
        setStandardYAxis(hAx2,[],6,yLim);
    end    
end

% resets the axis dimensions
resetAxesPos(hAx2,1,1);
resetObjPos(hAx2,'bottom',0.015,1);
resetObjPos(hAx2,'height',-0.030,1);

% resets the axis dimensions so that match the last sub-plot
dL = retObjDimPos(hAx2,1,1) - retObjDimPos(hAx,1,1);
resetObjPos(hAx2,'left',retObjDimPos(hAx,1,1));
if (pP.grpTime)    
    resetObjPos(hAx2,'width',retObjDimPos(hAx,3,1));    
else
    resetObjPos(hAx2,'width',dL,1);
    setGroupString(hAx2,pF,1:nGrp,lStr1,-90); 
end

% --- PRE/POST STIMULI SPEED MARKERS --- %
% -------------------------------------- %

% memory allocations
hPlot = zeros(nGrp,1);

% creates a new subplot
hAxM = subplot(2,2,3); hold on  
axis(hAxM,'off');

%
V1 = cellfun(@(x)(nanmean(x,1)),p.V1,'un',0);
V2 = cellfun(@(x)(nanmean(x,1)),p.V2,'un',0);
[V1,V2] = deal(V1(~cellfun(@isempty,V1(:))),V2(~cellfun(@isempty,V2(:))));
[V1,V2] = deal(cell2mat(V1(:)),cell2mat(V2(:)));

%
if (pP.plotPolar)
    D = sqrt(V1.^2 + V2.^2);
    Phi = atan2(V1,V2)*(180/pi);    
end

% loops through all the time groups plotting the speeds
for j = 1:nGrp
    if (sP.pT(j))
        [ii,jj] = deal(mod(j-1,ceil(nGrp/2))+1,floor((j-1)/ceil(nGrp/2))+1);
        if (~all(isnan(V1(:,j))))
            if (pP.plotPolar)
                [xPlt,yPlt] = deal(Phi(:,j),D(:,j));
            else
                [xPlt,yPlt] = deal(V1(:,j),V2(:,j));
            end
            
            % sets the plot values
            hPlot(j) = plot(hAxM,xPlt,yPlt,[c(jj),mm(ii)],'markersize',pP.mSize);                    
        end
    end
end    

% % updates the plot axis limits
if (pP.plotPolar)    
    set(hAxM,'xlim',[0 90],'box','on')    
else
    [V1T,V2T] = field2cell(plotD{1},{'V1','V2'});    
    if (pP.plotFixedY)     
        % determines the max mean V1 values
        V1Tmn = cellfun(@(y)(cellfun(@(x)(nanmean(x,1)),y(:),'un',0)),V1T,'un',0);
        V1Tmn = cellfun(@(y)(max(max(cell2mat(y(~cellfun(@isempty,y)))))),V1Tmn);
        
        % determines the max mean V2 values
        V2Tmn = cellfun(@(y)(cellfun(@(x)(nanmean(x,1)),y(:),'un',0)),V2T,'un',0);
        V2Tmn = cellfun(@(y)(max(max(cell2mat(y(~cellfun(@isempty,y)))))),V2Tmn);
        
        % sets the overall limit 
        YmxNw = detOverallLimit([detOverallLimit(V1Tmn) detOverallLimit(V2Tmn)]);
    else
        YmxNw = detOverallLimit([max(V1(:)),max(V2(:))]);
    end
    
    % sets the axis limits so that they are square
    set(hAxM,'xlim',[0 YmxNw],'ylim',[0 YmxNw],'box','on') 
end      

% updates the axis properties
formatPlotAxis(hAxM,pF,1); 
axis(hAxM,'normal');

% resets the axis location
hPos = get(hAxM,'Position');
resetObjPos(hAxM,'left',-0.15*hPos(3),1)
resetObjPos(hAxM,'width',0.20*hPos(3),1)                
    
% plots the gradient lines (if specified)
if (pP.addTrend)
    % calculates the font-size ratio
    fR = max(0.8,min(newSz(3:4)./regSz(3:4))*get(0,'ScreenPixelsPerInch')/72);
    fSzT = setMinFontSize(pF.Axis(1).Font.FontSize*fR,'text');
    
    % calculates the linear fits to the lines
    [M_mn,CI,R2] = calcLinearFit(V1,V2);    
    [xLim,yLim] = deal(get(hAxM,'xlim'),get(hAxM,'ylim'));

    % plots the values
    plot(hAxM,xLim,M_mn*xLim,'g','linewidth',lWid)
    plot(hAxM,xLim,CI(1)*xLim,'g--','linewidth',lWid)
    plot(hAxM,xLim,CI(2)*xLim,'g--','linewidth',lWid)        
    
    % sets the new text label
    nwText = [{sprintf('Gradient = %.2f %s %.2f',M_mn,'\pm',diff(CI)/2)};...
              {sprintf('(R^2 = %.2f)',R2)}];
    hText = text(0,0,nwText,'FontWeight',...
                 'bold','FontSize',fSzT,'HorizontalAlignment','center',...
                 'parent',hAxM,'tag','Other');
    
    [tExt,dX] = deal(get(hText,'Extent'),0.025*diff(xLim));
    set(hText,'position',[xLim(2)-(dX+tExt(3)/2),2*tExt(4)/3,0],'Units','Normalized')    
end

% creates the legend
if (nGrp > 1)
    [dY,isOut] = deal(0.02,sP.pT' & ~all(isnan(V1),1));
    if (sum(isOut) == 1)
        % only one valid time group
        hLg = legend(hPlot(isOut),lStr1(isOut),'location','NorthOutside');
        set(hLg,'box','off')
    else
        % more than one valid time group
        hLg = columnlegend(hPlot(isOut),2,lStr1(isOut),'location','NorthOutside');
    end

    % resets the axis locations
    [axPos,lgPos] = deal(get(hAxM,'position'),get(hLg,'position'));
        
    % resets the legend position
    lgPos(1) = axPos(1) + (axPos(3) - lgPos(3))/2;
    lgPos(2) = 1 - (dY + lgPos(4));
    set(hLg,'position',lgPos)
    
    % resets the axis height
    axPos(4) = (lgPos(2) + lgPos(4)/2) - (axPos(2) + 4*dY);
    set(hAxM,'position',axPos)
end

axis(hAxM,'on');

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)