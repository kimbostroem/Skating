function OutVect = calcKneeElbAxis(Vect,Angle,AltVect,TransitRng,isPlot)
% Handle elbow and knee joint vector when in near extension
% Use the AltVect with a smooth transition phase for the angular range TransitRng
% Vect         vector of jont axis as time series [3,nsamples]
% Angle        extension angle (180 deg is extended) [1,nsamples]
% AltVect      alternative vector (deg) [3,nsamples]
% TransitRng   transition range (deg) [min max]
% isPlot       (optional)
% version 1, 200615 (MdL)
% version 2, 210201 : resolve bug that occurred if isempty(AltVect)
% version 3, 210727 : resolve bug that prevented filling and the usage of the AltVect
%                     do not fill above the threshold
% version 4, 210802 : declare OutVect before possible return


% validity check
OutVect = nan(size(Vect));
if all(isnan(Vect),'all') && (isempty(AltVect) || all(isnan(AltVect),'all'))
    return;
end

Vect    = normLength(Vect);
AltVect = normLength(AltVect);

% Switch the regions of near extension with the AltVect across the TransitRange
Interp = Vect+(Angle-TransitRng(1))/(TransitRng(2)-TransitRng(1)).*(AltVect-Vect);
Small  = Angle<TransitRng(1);
Large  = Angle>TransitRng(2);
OutVect          = Interp;
OutVect(:,Small) = Vect(:,Small);
OutVect(:,Large) = AltVect(:,Large);
OutVect = normLength(OutVect);

% plotting if desired
if nargin >7 && isPlot
    figure; title('interpolate norm-vector in vicinity of full extension')
    subplot(2,2,1); hold on; plot(OutVect(1,:),'linewidth',1); plot(Vect(1,:));plot(AltVect(1,:));
    subplot(2,2,2); hold on; plot(OutVect(2,:),'linewidth',1); plot(Vect(2,:));plot(AltVect(2,:));
    subplot(2,2,3); hold on; plot(OutVect(3,:),'linewidth',1); plot(Vect(3,:));plot(AltVect(3,:));
    subplot(2,2,4); hold on; plot(Angle); yline(TransitRng(1)); yline(TransitRng(2)); title('angle');
end

end
