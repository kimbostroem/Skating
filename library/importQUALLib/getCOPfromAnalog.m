function [COPcorr,COP,Forces,IsLoadedAll,Type,Msgs,Flag] = getCOPfromAnalog(AnalogAll,Forces,FPNumbers,Locations,Freq, ...
    FiltFreq,PowerLineFreq,Fill,COPraw,Type,DoCrop)
%% Smooth and correct the COP, depending on the type of force plate, using the original analog data if present
% corrected and filtered COP for Kistler 9287CA force plates
% since we are dealing with bit-noise, the noise depends on amplification.
% therefore it is best to apply the load threshold on the analog z-channels
% 
% SYNTAX
% [COPcorr,COP,Forces,IsLoadedAll,Type,Msgs,Flag] = getCOPfromAnalog(AnalogAll,Forces,FPNumbers,Locations,Freq,FiltFreq,PowerLineFreq,Fill,COPraw,Type,DoCrop)
% 
% INPUT
%   AnalogAll   (double Plates x 8-Chan x NSamples) timeseries of the eight analog channels
%   Forces      (double Plates x 3d x NSamples) timeseries of the 3D forces
%   FPno        (double Plates) list of Force plate numbers to be analysed (default: all eight plates)
%   Locations   (double Plates x 4-corners x 3d) locations of the corners of the force plates [m, meter]
%   Freq        (double) measurement frequency [Hz]
% Optional:
%   FiltFreq    (double) [default 20 Hz; 0=off] low-pass filter depending on Filter:
%   PowerFreq   (double) [default 50 Hz; 0=off] powerline hum removal
%   Fill        (double) [default NaN] COP value for non-loaded phases (0=force plate center)
%   COPraw      (double) If Analog is not Kistler-type, the COP is not re-calculated, only smoothed.
%   Type        (cell char) for Kistler FP, the COP is calculated directly from analog signals
%   DoCrop      (logical) default true: crop the forces when outside the force plate edges
%            
% OUTPUT
%     COPcorr     (Plates x 3d x NSamples double) COP corrected for non-linearity of force plates 
%     COP         (Plates x 3d x NSamples double) COP not corrected for non-linearity of force plates
%     Forces      (Plates x 3d x NSamples double) forces with removed Powerline hum and gap-filled
%     IsLoadedAll (Plates x NSamples logical) which of te plates is loaded when
%     Type        (char) type of force plate 
%     Flag        (double) Exit status of the function (1 = success, 0 = warning, -1 = error)
%     Msgs        (char) notifications, warning and error messages
%
% EXAMPLES
% - Use with defaults:
% COP = getCOPfromAnalog(AnalogAll,Forces,[],Locations,Freq);
% 
% - Using only a selection of Force Plates, and no filtering:
% COP = getCOPfromAnalog(AnalogAll,Forces,[1 5],Locations,Freq,0,0);
% 
% - In case AMTI plates are also used:
% COP = getCOPfromAnalog(AnalogAll,Forces,[],Locations,Freq,[],[],[],COPQTMall);
%
% Local functions: kistlerCOPFromAnalog, correctCOP6d
% 
% See also: computeLoadedPhases
% 
% (c) 2019 by WWU Muenster
% Marc de Lussanet
% version 230204 (MdL) fixed condition for reconstruction of kistler COP
% version 230320 (MdL) fixed erroneous restructuring if there is just one force plate
% version 230516 (MdL) fixed wrong loading thresholds for Kistler
% version 230803 (MdL) flag for cropping the forces to the force plate edges; DebugPlot; cleaning of header

%% Defaults
narginchk(5,11)
if nargin < 6    || isempty(FiltFreq),      FiltFreq     = 20;   end
if nargin < 7    || isempty(PowerLineFreq), PowerLineFreq= 50;   end
if nargin < 8    || isempty(Fill),          Fill         =  0;   end
if nargin < 9    || isempty(COPraw),        COPraw       = [];   end
if nargin < 10   || isempty(Type),          Type         = {''}; end
if nargin < 11   || isempty(DoCrop),        DoCrop       = false; end
Msgs = {};
Flag = 1;
IsPLOT = false;

% The filter frequency for the COP must be not too low to deal with the transients in the Z
% direction (e.g. jumps)
FiltFreq_COP = max([FiltFreq 40]); 

% constants for determination of loading
MinimalLoadDuration = 0.1 * Freq; % defailt 0.1 s : very short loaded periods are not plausible
LoadThresholds      = [20 10];    % N: rough and fine estimate for load detection

%% Constants & definitions
X=1; Y=2; Z=3;
NPlates     = size(Forces,1);    % Number of plates
NSamples    = size(Forces,3);
COP         = Fill * ones(NPlates,Z,NSamples); % COP computed from analog signals
COPcorr     = COP; % COP corrected for Kistler plate deformation
IsLoadedAll = false(NPlates,NSamples);
% read only force plates that are requested and which are present in the data
if isempty(FPNumbers)
    ForcePlates = 1:NPlates;
else
    FPNumbers(FPNumbers>NPlates) = [];
    ForcePlates = FPNumbers;
end

%% loop through the force plates
for iFP = ForcePlates
    % init
    Force    = squeeze(Forces(   iFP,:,:));
    Location = squeeze(Locations(iFP,:,:));
    
    %% remove power hum (filtering should be on analog rather than COP data!)
    if PowerLineFreq && not(any(contains(Type,'AMTI'))) % AMTI force plates do not have power line hum
        Force = periodicMedianFilter(Force,   round(Freq/PowerLineFreq));
    end
    
    %% find non-loaded phases and set to nan
    % boolean array for each sample if the current plate is loaded
    IsLoaded = computeLoadedPhases(Force(Z,:), LoadThresholds,MinimalLoadDuration,iFP);
    IsLoadedAll(iFP,:) = IsLoaded;

    %% Compute the COP
    if any(IsLoaded) % (only if the plate was loaded at all)
        %% low pass filtering if desired
        % fill non-loaded periods with NaN values
        Force(  :,~IsLoaded) = 0;
        if FiltFreq == 0
            % do not filter
        else % 'Butter'
            ForceForCOP = gapfilth(FiltFreq_COP,Freq,2,Force,   'l',2);
            Force   = gapfilth(FiltFreq,Freq,2,Force,   'l',2);
        end
        
        %% Calculate a smooth COP
        DoGetCOPFromAnalog = any(contains(Type,'Kistler')) && ~isempty(AnalogAll);
        if DoGetCOPFromAnalog
            Analog   = squeeze(AnalogAll(iFP,:,:));
            [COPnoCorr,COPno] = kistlerCOPFromAnalog(Analog,ForceForCOP,IsLoaded,Location,iFP,Freq,PowerLineFreq,FiltFreq_COP);

        elseif ~isempty(COPraw) % only if COPraw is provided 
            COPno = squeeze(COPraw(iFP,X:Y,:));
            COPno(:,~IsLoaded)=NaN;
            if FiltFreq == 0
                % do not filter
            else % 'Butter'
                COPno = gapfilth(FiltFreq_COP,Freq,2,COPno,'l',2);
            end
            COPnoCorr = COPno;
        else
            error('the COPraw must be provided if the focre is not provided from Qualisys-Kistler measurements')
        end

        % Set non-loaded phases to Fill if wanted (Fill==0 : force plate center)
        COPno(    :,~IsLoaded) = Fill;
        COPnoCorr(:,~IsLoaded) = Fill;

        % copy the results of the current plate
        COP(iFP,X:Y,:) = COPno;
        COPcorr(iFP,X:Y,:) = COPnoCorr;
        Forces(iFP,:,:) = Force;
    else
        % Set non-loaded phases to Fill if wanted (Fill==0 : force plate center)
        COP(iFP,X:Y,:) = Fill;
        COPcorr(iFP,X:Y,:) = Fill;
        Forces(iFP,:,:) = 0;
    end

    %% shift to Forceplate location
    if any(contains(Type,'Kistler'))
        LcX = mean(Location(:,X)); % plate center in X
        LcY = mean(Location(:,Y)); % plate centre in Y
        COP( iFP,X,:) = COP( iFP,X,:) + LcX;
        COP( iFP,Y,:) = COP( iFP,Y,:) + LcY;
        COPcorr(iFP,X,:) = COPcorr(iFP,X,:) + LcX;
        COPcorr(iFP,Y,:) = COPcorr(iFP,Y,:) + LcY;
    end

    %% if just one plate is requested, return the data just for that plate
    if length(FPNumbers)==1
        COPcorr = squeeze(COPcorr(iFP,:,:));
        COP     = squeeze(COP(iFP,:,:));
    end
end

%% crop to the edges of the force plate if desired
if DoCrop
    COP     = cropForcesToForcePlates(COP,    [],Locations,Fill);
    COPcorr = cropForcesToForcePlates(COPcorr,[],Locations,Fill);
end

% Debug plot if desired
if IsPLOT
    figure;hold on;
    title('getCOPfromAnalog : the calculated COP for all plates')
    for iFP=1:NPlates
        COPi=squeeze(COPcorr(iFP,:,:));
        COPi(:,~IsLoadedAll(iFP,:)) = nan;
        plot3(COPi(1,:),COPi(2,:),COPi(3,:),'.-');
        Loc = squeeze(Locations(iFP,:,:));
        Loc = [Loc; Loc(1,:)]; 
        plot3(Loc(:,1),Loc(:,2),Loc(:,3),'linewidth',2)
    end
end
end

%% ========================================================================

function [COPcorr,COP] = kistlerCOPFromAnalog(Analog,Force,IsLoaded,Location,FPNo,Freq,PowerLineFreq,FiltFreq)
%% Calculate a smooth COP

% Kistler parameters
KistlerDepth  =  0.053; % [m]
KistlerRadius = [0.350 0.210];
Calibration   = [ ...
    38.472,38.258,38.118,38.042, 19.509,19.494,19.587,19.370; ...
    38.234,38.256,38.236,37.954, 19.317,19.155,19.140,19.387; ...
    38.199,38.131,38.067,38.273, 19.523,19.477,19.469,19.411; ...
    37.817,37.733,37.803,38.017, 19.326,19.369,19.350,19.295; ...
    38.493,38.147,38.304,38.036, 19.407,19.516,19.419,19.318; ...
    38.695,38.465,38.406,38.514, 19.303,19.360,19.391,19.466; ...
    38.423,38.577,38.532,38.418, 19.332,19.390,19.390,19.448; ...
    38.403,38.517,38.010,38.210, 19.579,19.796,19.529,19.576];

% 
X=1;Y=2;Z=3;
ZChannels = 5:8;
AnalogZ = squeeze(Analog(ZChannels,:)) ./ Calibration(FPNo,ZChannels)';  % Z direction of analog signals
if size(AnalogZ,2)==1
    AnalogZ = AnalogZ'; % bring to the form 1xN
end
AnalogZ = periodicMedianFilter(AnalogZ, round(Freq/PowerLineFreq));

%% low pass filtering if desired
% fill non-loaded periods with NaN values
AnalogZ(:,~IsLoaded) = NaN;
if FiltFreq == 0
    % do not filter
else % Butterworth filter
    AnalogZ = gapfilth(FiltFreq,Freq,2,AnalogZ, 'l',2); %#ok<*UNRCH>
end

%% calculate COP from analog channels (only vertical)
SumAna = sum(AnalogZ);
% this is the normalised location of the COP
COP = [(-AnalogZ(1,:)-AnalogZ(2,:)+AnalogZ(3,:)+AnalogZ(4,:)) ./ SumAna ; ...
    (-AnalogZ(1,:)+AnalogZ(2,:)+AnalogZ(3,:)-AnalogZ(4,:)) ./ SumAna ];
% convert to meters
SignX = sign(Location(3,X)-Location(1,X));    % some plates are turned 180 deg
SignY = sign(Location(2,Y)-Location(1,Y));    % some plates are turned 180 deg
COP(X,:) = COP(X,:) .* KistlerRadius(X) * SignX; % sensors  are 10 cm from x-edge
COP(Y,:) = COP(Y,:) .* KistlerRadius(Y) * SignY; % sensors  are  9 cm from y-edge

%% correct for oblique force direction (Dominik.Jenni@kistler.com: 52 mm, niet 53)
% We use the force, not the Analog, because Z can have a different gain than XY
COP(X,:) = COP(X,:) + KistlerDepth * Force(X,:) ./ Force(Z,:);
COP(Y,:) = COP(Y,:) + KistlerDepth * Force(Y,:) ./ Force(Z,:);

%% The COP must be inside the plate surface
Size = (max(Location,[],1)-min(Location,[],1))/2;
Out = abs(COP(X,:))>Size(X) | abs(COP(Y,:))>Size(Y);
COP(:,Out) = NaN;

%% correct errors by variable polynomial (version 4 by Tobias: 16.7.2019)
COPcorr = correctCOP6d(COP);  % fits slightly better
end

%% ========================================================================

function [COPcorr] = correctCOP6d(COP)
%%input: COP of force plate Kistler 9287CA
%% values should be between x -0.45:0.45 and y -0.3:0.3
%% Error correction using curve fitting tool polynomial 4th order
%% Fitting parameters for Plate no 1 (SN4327853), measurement of June 2019

X = 1; Y = 2;
x = COP(X,:); y = COP(Y,:);

% X error Coefficients: symmetry demands only odd powers for x and even ones for y
a0 = -15210;
a1 =   3051;
a2 =   -159.6;
a3 =      1.304;
b0 =   5980;
b1 =  -1014;
b2 =     44.78;
b3 =     -0.7626;
c0 =   -470.6;
c1 =     67.98;
c2 =     -2.261;
c3 =      0.04445;
Cx= (a0*y.^6 + a1*y.^4 + a2*y.^2 + a3) .*x.^5 + (b0*y.^6 + b1*y.^4 + b2*y.^2 + b3) .*x.^3 + (c0*y.^6 + c1*y.^4 + c2*y.^2 + c3) .*x;

% Y error Coefficients: symmetry demands only odd powers for y and even ones for x
a0 = 21440;
a1 = -6481;
a2 =   419.3;
a3 =     5.212;
b0 = -2041;
b1 =   723.9;
b2 =   -63.39;
b3 =    -0.7349;
c0 =    23.6;
c1 =   -16.12;
c2 =     2.487;
c3 =    -0.05381;
Cy= (a0*x.^6 + a1*x.^4 + a2*x.^2 + a3) .*y.^5 + (b0*x.^6 + b1*x.^4 + b2*x.^2 + b3) .*y.^3 + (c0*x.^6 + c1*x.^4 + c2*x.^2 + c3) .*y;

COPcorr      = COP;
COPcorr(X,:) = COP(X,:) - Cx;
COPcorr(Y,:) = COP(Y,:) - Cy;

end



