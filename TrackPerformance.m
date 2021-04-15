% --- diagnostic function that illustrates the efficacy of the tracking 
%     process
function pData = TrackPerformance(snTot)

% if no plot data struct is given, then initialise a new one
pData = initPlotDataStruct(mfilename,@calcFunc,@plotFunc,@outFunc);    

% sets the function name/type
pData.Name = 'Overall Tracking Performance';
pData.Type = 'Indiv';
pData.fType = [2 1 1 1];

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
    [pData.hasTime,pData.hasSR] = deal(true,true);
    pData.sP = initSpecialPara(snTotL,pData,[],1); 
end
    
% ----------------------------------------------------------------------- %
% ---                 PARAMETER STRUCT SETUP FUNCTIONS                --- %
% ----------------------------------------------------------------------- %

% --- initialises the calculation parameter function --- %
function cP = initCalcPara(snTot)

% initialises the parameter struct
nPara = 0;
cP = setParaFields(nPara);

% --- initialises the calculation parameter function --- %
function pP = initPlotPara(snTot)

% initialises the parameter struct
nPara = ~detIfShortExpt(field2cell(snTot,'T'));                       
pP = setParaFields(nPara);

% sets the tab list names
a = {'1 - General'};

% sets the plot parameter fields into the data struct
if (nPara > 0)
    pP(1) = setParaFields(a{1},'Boolean',0,'isZeitG','Use Zeitgeiber Time Format');
end

% --- initialises the plot formatting data struct --- %
function pF = initPlotFormat(snTot)

% memory allocation
nApp = length(snTot.appPara.ok); 
pF = setFormatFields(1);

% initialises the font structs
pF.Title = setFormatFields([],'',nApp);
pF.xLabel = setFormatFields([],' ',1);
pF.yLabel = setFormatFields([],'Fly Index',1);
pF.Axis = setFormatFields([],[]);

% sets the apparatus names as the titles
for i = 1:nApp
    pF.Title(i).String = snTot.appPara.Name{i};
end

% --- initialises the output data parameter struct --- %
function oP = initOutputPara(snTot)

% initialisations
oP = [];

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

% other initialisations
flyok = snTot.appPara.flyok;
nApp = length(flyok);
T = cell2mat(snTot.T);

% parameters
[dTol,nFrm] = deal(20,length(T));

% memory allocation
[plotD,ok] = deal(repmat(struct('T',T,'Y',[]),nApp,1),false);

% creates the waitbar figure
h = ProgBar({'Overall Progress'},'Performance Tracking Calculations');

% ---------------------------- %
% --- FLY LOCATION SETTING --- %
% ---------------------------- %

% array dimensioning
flyok = snTot.appPara.flyok;

% sets the position values for all the apparatus
for i = 1:nApp
    % updates the waitbar figure
    wStrNw = sprintf(['Setting Up Performance Image Arrays ',...
                      '(Region %i of %i)'],i,nApp);
    if h.Update(1,wStrNw,(i+1)/(2+nApp))
        % if the user cancelled, then exit the function
        [plotD,ok] = deal([],false);
        return
    end       
    
    % sets the y-coordinates (if they exist)
    if ~isempty(snTot.Px)
        Px = snTot.Px{i}';
    else
        Px = zeros(length(flyok{i}),nFrm);
    end
    
    % sets the y-coordinates (if they exist)
    if ~isempty(snTot.Py)
        Py = snTot.Px{i}';
    else
        Py = zeros(length(flyok{i}),nFrm);
    end   
    
    % calculates the inter-frame displacements over all flies and frames
    [dX,dY] = deal(diff([Px(:,1),Px],[],2),diff([Py(:,1),Py],[],2));
    D = sqrt(dX.^2 + dY.^2);
    
    % array initialisation
    Y = zeros([size(D),3]);
    Y(:,:,2) = 255;
    
    % sets the non-NaN indices        
    for j = 1:size(D,1)
        if ~flyok{i}(j)
            % flags that the fly has been rejected
            Y(j,:,:) = 125;
        else
            % otherwise, determine the rejected/large distance frames
            indR = min(length(T),roundP(T(isnan(D(j,:)))+1));
            indB = min(length(T),roundP(T(D(j,:) > dTol)+1));

            % set the frames in the overall array
            [Y(j,indR,1),Y(j,indR,2)] = deal(255,0);
            [Y(j,indB,3),Y(j,indB,2)] = deal(255,0);
        end
    end
    
    % sets the final array
    plotD(i).Y = Y;
end

% closes the waitbar figure
if ~h.Update(1,'Performance Image Array Setup Complete!',1)
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

% memory allocation and other dimensioning
p = plotD{1}(sP.pInd);
if (isempty(p.Y)); return; end

% sets the plot indices
ind0 = find(p.T >= sP.xLim(1),1,'first');
indF = find(p.T <= sP.xLim(2),1,'last');

% other parameters
nFly = size(p.Y,1);

% retrieves the panel object handle
hP = get(gca,'Parent');

%
if (isfield(pP,'isZeitG'))
    isZeitG = pP.isZeitG;
else
    isZeitG = false;
end

% ---------------------------------------- %
% --- FORMATTING STRUCT INTIALISATIONS --- %
% ---------------------------------------- %

% reformats the font data struct
pF = retFormatStruct(pF,1);

% resets the axis font size if greater than threshold
pF.Axis(1).Font.FontSize = pF.Axis(1).Font.FontSize - 4*(nFly > 20);

% case is using the absolute time axis
if (isZeitG)            
    % case is using Zeitgeiber time
    pF.xLabel.String = 'Zeitgeiber Time';
else
    % if plotting day/night, then set absolute time
    pF.xLabel.String = 'Time'; 
end

% ----------------------- %
% --- FIGURE CREATION --- %
% ----------------------- %

% sets the axis handle
hAx = createSubPlotAxes(hP);
hold(hAx,'on'); axis(hAx,'on')

% plots the image
imagesc(p.Y/255);

% case is using the absolute time axis
if (isZeitG)            
    % case is using Zeitgeiber time
    pF.xLabel(1).String = 'Time (Zeitgeiber)';
    setZeitGTimeAxis(hAx,p.T,snTot);        
else
    % if plotting day/night, then set absolute time
    pF.xLabel(1).String = 'Time';
    setAbsTimeAxis(hAx,p.T,snTot);         
end

% sets the axis properties
% xLim = p.T([ind0,indF]);
xLim = [ind0,indF];
set(hAx,'xlim',xLim,'ylim',[1 nFly]+0.5*[-1 1],'ytick',1:nFly,...
        'yticklabel',num2cell(1:nFly),'xgrid','on')

% resets the axis positions
axis(hAx,'ij')

%
formatPlotAxis(hAx,pF,sP.pInd);
resetAxesPos(hAx,1,1); 

% if the plots do not exist, then create them for each apparatus
for i = 1:nFly-1    
    % adds the seperator markers
    plot(xLim,(i+0.5)*[1 1],'k','linewidth',1)
end

% ----------------------------------------------------------------------- %
% ---                         OUTPUT FUNCTION                         --- %
% ----------------------------------------------------------------------- %

% --- function for calculating the data for the analysis figure --- %
function [pData,plotD] = outFunc(pData,plotD,snTot)