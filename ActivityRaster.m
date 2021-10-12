% --- plots the rasterplot of the individual fly activity over the duration 
%     the experiment (short experiment only)
function pData = ActivityRaster(snTot)

% initialises the plot data struct
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Activity Rasterplot (Short)';
pData.Type = {'Indiv','Multi'};
pData.fType = [1 1 2 1];
pData.rI = initFuncReqInfo(pData);

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
    pData.hasSP = true;
    pData.sP = initSpecialPara(snTotL,pData); 
end
    
% ----------------------------------------------------------------------- %
% ---                 PARAMETER STRUCT SETUP FUNCTIONS                --- %
% ----------------------------------------------------------------------- %

% --- sets the function required information struct
function rI = initFuncReqInfo(pData)

% memory allocation
rI = struct('Scope',[],'Dur',[],'Shape',[],...
            'Stim',[],'Spec',[],'SpecFcn',[]);

% sets the struct fields
rI.Scope = setFuncScopeString(pData.Type);
rI.Dur = 'Short';
rI.Shape = 'None';
rI.Stim = 'None';
rI.Spec = 'None';

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 4;                       
cP = setParaFields(nPara);

% sets the parameter fields
cP(1) = setParaFields([],'Number',3,'vAct','Activity Threshold (mm/s)',[0.1 10 false]);
cP(2) = setParaFields([],'Boolean',1,'useAll','Analyse Entire Experiment');
cP(3) = setParaFields([],'Number',0,'T0','Start Time (min)',[0 inf true],{2,1});
cP(4) = setParaFields([],'Number',30,'Tdur','Analysis Duration (min)',[10 inf true],{2,1});

% sets the tool-tip strings
cP(1).TTstr = 'The inter-frame velocity threshold used to indicate movement';
cP(2).TTstr = 'Determines whether or not to use the entire experiment for the analysis';
cP(3).TTstr = 'The time into the experiment where the analysis begins';
cP(4).TTstr = 'The duration period of the analysis';

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% retrieves the stimuli marker time-stamps
Ts = arrayfun(@(x)(getMotorFiringTimes(x.stimP)),snTot,'un',0);
hasStim = ~any(cellfun(@isempty,Ts));

% initialises the parameter struct
nPara = 4 + hasStim;                      
pP = setParaFields(nPara);

% tab list headings
a = {'1 - General','2 - Raster Plot'};

% sets the parameter fields
pP(1) = setParaFields(a{2},'Number',10,'nHght','Raster Plot Bar Height',[1,20,true]);
pP(2) = setParaFields(a{2},'Number',2,'nGap','Raster Plot Gap Size',[0,10,true]);
pP(3) = setParaFields(a{1},'Boolean',false,'showNum','Show Fly Index Numbering');
pP(4) = setParaFields(a{1},'Boolean',false,'setEqual','Set Equal Index Counts');

% only set this field if there are stimuli events present in expt
if hasStim
    pP(5) = setParaFields(a{1},'Boolean',true,'showStim','Show Stimuli Markers');
end

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% determines the apparatus count
nApp = length(snTot.iMov.ok);   
pF = setFormatFields(nApp);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],'',1);
pF.yLabel = setFormatFields([],'',1);
pF.Axis = setFormatFields([],[]);

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Title(i).String = snTot.iMov.pInfo.gName{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = setupOutputParaStruct(snTot,false,true,false);

% sets the independent variable fields
oP = addXVarField(oP,'Time','TT','Time');

% sets the dependent output variables
oP = addYVarField(oP,'Activity','I',[],5,{'TT'},1);
                         
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
Tf = cellfun(@(x)(x(end)),T)/60;

% checks to see if the solution struct has the sub-region data struct
ok = checkFuncPara({'AnalysisTimeCheck'},cP,Tf);
if (~ok); plotD = []; return; end

% sets the movement calculation type
cP.movType = 'Absolute Speed';

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% array dimensions
nExp = length(snTot);

% memory allocation
nApp = length(snTot(1).iMov.ok);
plotD = initPlotValueStruct(snTot,pData,cP,'T',[],'TT',[],'I',[],'IT',[]);                        
[I,TT] = deal(cell(nApp,nExp),cell(1,nExp));

% creates the waitbar figure
wStr = {'Overall Progress'};
h = ProgBar(wStr,'Activity Calculations');

% ------------------------------------------------------- %
% --- INTER-FRAME DISTANCE CALCULATION & THRESHOLDING --- %
% ------------------------------------------------------- %
    
% loops through each of the experiments calculating the velocity values
for i = 1:nExp 
    % updates the waitbar figure (if more than one solution file)
    wStrNw = sprintf('%s (Experiment %i of %i)',wStr{1},i,nExp);
    if h.Update(1,wStrNw,i/(1+nExp))
        [plotD,ok] = deal([],false);
        return
    end

    % calculates the video frame rate and experiment apparatus indices
    FPS = snTot(i).sgP.fRate/snTot(i).sgP.sRate;
    iApp = find(~cellfun(@isempty,snTot(i).iMov.flyok));
    
    % sets the relevant time points and apparatus indices for this expt
    if cP.useAll
        % uses all the time points
        cP.T0 = 0;
        ii = 1:length(T{i});
    else
        % use only the points from the start to the duration end
        ii = (T{i} >= 60*cP.T0) & (T{i} <= 60*(cP.T0 + cP.Tdur));
    end      
    
    % sets the relevant time points and apparatus indices for this expt
    % (after the start time and before the end of the analysis duration) 
    Tnw = T{i}(ii)-60*cP.T0;     
    isMove = calcFlyMove(snTot(i),Tnw,ii,iApp,cP.vAct);      
    
    % calculates the time mid point of each frame, and removes all the
    % frames where the inter-frame time difference is large (ie, the
    % change-over in video recordings)
    [Tmid,isOK] = deal((Tnw(1:end-1) + Tnw(2:end))/120,diff(Tnw) < 2/FPS);
    
    % sets the 
    TT{i} = Tmid(isOK);
    for j = 1:length(isMove)
        fok = snTot(i).iMov.flyok{iApp(j)};
        plotD(iApp(j)).I(1,fok,i) = num2cell(isMove{j}(isOK,:),1);
        plotD(iApp(j)).TT{i} = TT{i}*60;
    end
end
    
% sets raster plot array for each apparatus
for i = 1:nApp
    % sets the inactivity values for the current group
    Inw = plotD(i).I(:)';
    Inw = Inw(~cellfun(@isempty,Inw));
            
    % sets all of the experiments into a single array
    plotD(i).IT = false(max(cellfun(@length,Inw)),length(Inw));
    for j = 1:length(Inw)
        plotD(i).IT(1:length(Inw{j}),j) = Inw{j};
    end
    
    % sets the time vector
    plotD(i).T = TT{argMax(cellfun(@length,TT))}*60;
end
    
% closes the waitbar figure
if ~h.Update(1,'Activity Calculations Complete!',1)
    if nargin < 5 
        h.closeProgBar(); 
        pause(0.05); 
    end
end

% ----------------------------------------------------------------------- %
% ---                        PLOTTING FUNCTION                        --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,varargout] = plotFunc(snTot,pData,plotD,ind)

% retrieves the plotting paraeter struct
cP = retParaStruct(pData.cP);
pP = retParaStruct(pData.pP);
sP = retParaStruct(pData.sP);
pF = pData.pF;

% creates a loadbar
h = ProgressLoadbar('Creating Activty Rasterplot...');

% ------------------------------------------- %
% --- INITIALISATIONS & MEMORY ALLOCATION --- %
% ------------------------------------------- %

% sets the plotting indices and subplot indices
[ind,m,n] = deal(find(sP.Sub.isPlot),sP.Sub.nRow,sP.Sub.nCol);
nApp = length(ind); if (nApp == 0); return; end
p = plotD{1}(ind);

% retrieves the stimuli time stamps
if isfield(snTot,'Ts')
    Ts = cellfun(@(x)(cell2mat(x)/60),field2cell(snTot,'Ts'),'un',0);
else
    Ts = cell(size(snTot));
end

% converts the solution time arrays into single vectors
if (cP.useAll)
    T = cellfun(@(x)(cell2mat(x)),field2cell(snTot,'T'),'un',0);
    Tf = cellfun(@(x)(x(end)),T)/60;
	[cP.T0,cP.Tdur] = deal(0,max(Tf));
end

% sets the other indices/parameters
nFrm = cellfun(@(x)(size(x,1)),field2cell(p,'IT'));
[nFly,nFlySz] = deal(cellfun(@(x)(size(x,2)),field2cell(p,'IT')));
if (pP.setEqual)
    nFlySz = max(nFly)*ones(size(nFly));    
end

% other parameters
pGray = 0.8;

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% titles are empty, so reset them to the apparatus names
for i = 1:length(pF.Title)
    pF.Title(i).String = snTot(1).iMov.pInfo.gName{i}; 
end

% sets the output data
pData.pF = pF;

% retrieves the formatting struct
if (isempty(m)); szMx = 1; else szMx = max([m n]); end
pF = retFormatStruct(pData.pF,szMx);

% sets the x/y labels
[pF.xLabel.String,pF.yLabel.String] = deal('Frame Index','Fly Index');
pF.Title = pF.Title(ind);

% calculates the y-axis tick markers
yy = (pP.nHght+pP.nGap); 

% determines the frame index when the stimuli events took place
ii = ~cellfun(@isempty,Ts);
if (any(ii))
    TsMn = groupStimTimes(cell2mat(Ts));
    TsMn = TsMn((TsMn - cP.T0) < cP.Tdur) - cP.T0;        
    iStim = cellfun(@(x)(argMin(abs(p(1).T-60*x))),num2cell(TsMn));
end

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% if nothing selected, then exit the function
hAx = cell(1,nApp);

% loops through all the subplots 
for j = 1:nApp
    % sets the actual plot index and creates the subplot region
    [hAx{j},k] = deal(subplot(m,n,j),ind(j));  
    [pF.xLabel.ind,pF.yLabel.ind] = deal(NaN); 
    
    % sets up the plot region
    I = ones(nFlySz(j)*pP.nHght + (nFlySz(j)+1)*pP.nGap,nFrm(j));
    for i = 1:nFly(j)
        ii = (i-1)*(pP.nHght + pP.nGap) + ((1:pP.nHght) + pP.nGap);
        I(ii,:) = repmat(1-p(j).IT(:,i)',pP.nHght,1);
    end
    I((ii(end)+pP.nGap+1):end,:) = pGray;
    
    % creates the new image and sets the image properties
    imagesc(I); axis normal; axis ij; hold on
    colormap('gray');    
    
    % plots the SEM filled regions (if set)
    if isfield(pP,'showStim')
        if pP.showStim
            for i = 1:length(iStim)
                plot(iStim(i)*[1 1],[0 size(I,1)]+0.5,'r','linewidth',3)
            end
        end
    end
    
    % sets the y-axis label strings
    yTick = (0:yy:(nFly(j)-1)*yy) + (pP.nHght/2+pP.nGap+1);
    if (pP.showNum)
        % sets the tick strings to the fly indices
        yTickLbl = num2cell(1:nFly(j))';
    else
        % sets the tick strings to be empty
        yTickLbl = repmat({''},nFly(j),1);
    end    
        
    % sets the x/y axis limits    
    set(hAx{j},'TickLength',[0 0],'ytick',yTick,'yticklabel',yTickLbl,...
               'xlim',[0 nFrm(j)+1],'UserData',j)
    formatPlotAxis(hAx{j},pF,j);        
end

% sets the non-aligned x/y labels
formatMultiXYLabels(hAx,pF,[m,n]);   

% % resets the object positions
% resetAxesPos(hAx,m,n);

% --- PLOT AXES REFORMATTING --- %
% ------------------------------ %

% simple error fix - plots the x-axis as it doesn't appear to be turning up
% automatically for the MacOS version?
if (ismac)
    for i = 1:length(hAx)    
        plot(get(hAx{i},'xlim'),(min(get(hAx{i},'ylim'))+eps)*[1 1],'k','linewidth',1.5) 
    end
end

% creates a loadbar
pause(0.05);
try; delete(h); end
    
% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)

% ----------------------------------------------------------------------- %
% ---                         OTHER FUNCTIONS                         --- %
% ----------------------------------------------------------------------- %

% --- determines the grouped stimuli times
function TsMn = groupStimTimes(Ts)

% tolerance
[tol,k] = deal(0.5,2);

% determines if the 
dTs = abs(Ts - mean(Ts));
if (any(dTs > tol))
    % if so, use k-means to cluster the data
    while (1)
        % calculates the k-means clustering for the current k-value
        ii = kmeans(Ts, k);        
        
        % calculates the mean differences within the groups
        [dTs,TsMn] = deal(zeros(k,1));
        for i = 1:k
            jj = (ii == i);
            TsMn(i) = mean(Ts(jj));
            dTs(i) = max(abs(Ts(jj) - TsMn(i)));
        end
        
        % determines if the differences are greater than tolerance
        if all(dTs < tol)
            % if not, then exit the loop
            break
        else
            % otherwise, continue the search
            k = k + 1;
        end
    end
else
    % otherwise, return the mean value
    TsMn = mean(Ts);
end
