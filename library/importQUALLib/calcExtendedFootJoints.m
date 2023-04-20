function Segments = calcExtendedFootJoints(Segments, Markers)
%% Create the five-segmented foot according to the Ghent foot model
% Calculate the Ghent kinematic foot model
% A Clinically Applicable Six-Segmented Foot Model (Sophie De Mits,1 Veerle Segers,2 Jim Woodburn,3
% Dirk Elewaut,4 Dirk De Clercq,2 Philip Roosen) JOURNAL OF ORTHOPAEDIC RESEARCH APRIL 2012
% 
% Usage:
%     SimInMN = importQUAL(FPath);
%     [~,EstBody] = importQUAL(FPath,~,DoEstLengths,RefBody);
%
% input-parameters:
%    Segments   (struct) Coordinate systems of the segments
%     Markers   (struct of 3 x N) time series for each marker that is present
%
% output-parameters:
%    Segments   (struct) Coordinate systems of the segments
% 
% See also: importQUAL, calcMainCS
%
% (c) Movement Science, WWU Muenster
% Authors: Rosemary Dubbeldam
% 230208 (MdL) outsourced isempty(Segments); updated parmeter list

%  Add Shank and FootV3d segments RDub
rAnkleVector= normLength(Markers.right_med_malleolus_marker-Markers.right_lat_malleolus_marker);
lAnkleVector= normLength(Markers.left_lat_malleolus_marker -Markers.left_med_malleolus_marker);
rForefootVector= normLength(Markers.right_DM1_marker-Markers.right_DM5_marker);
lForefootVector= normLength(Markers.left_DM5_marker -Markers.left_DM1_marker);
rKneeJoint  = squeeze(Segments.rTibia(4,:,:));
lKneeJoint  = squeeze(Segments.lTibia(4,:,:));
rAnkleJoint = squeeze(Segments.rFoot(4,:,:));
lAnkleJoint = squeeze(Segments.lFoot(4,:,:));
 
% Add V3D Shank coordinate system RDub
Segments.rShank    = createCS(rKneeJoint, rAnkleJoint-rKneeJoint, -rAnkleVector,'zxy');
Segments.lShank    = createCS(lKneeJoint, lAnkleJoint-lKneeJoint, -lAnkleVector,'zxy');
%
Segments.rFootV3D  = createCS(rAnkleJoint,(Markers.right_DM1_marker+Markers.right_DM5_marker)/2 - rAnkleJoint, -rForefootVector,'zxy');
Segments.lFootV3D  = createCS(lAnkleJoint,(Markers.left_DM1_marker+Markers.left_DM5_marker)/2   - lAnkleJoint,  -lForefootVector,'zxy');

% right foot segments RDub
rhfjc = (Markers.right_calcaneus_marker + Markers.right_calcaneus1_marker )/2; % right hindfoot joint centre
rmfjc = (Markers.right_NAV_marker + Markers.right_CUB_marker )/2;  % right midfoot joint centre
rmfdc = (Markers.right_PM1_marker + Markers.right_PM5_marker )/2;  % right midfoot distal centre
% virtual/additional points
RVPM2   = ( 2 * Markers.right_PM1_marker + Markers.right_PM5_marker )/3; % right virtual PM2 marker
baselineRDM = ( Markers.right_DM5_marker - Markers.right_DM1_marker );   % right distal forefoot baseline
baselineRDM = normLength(baselineRDM);

% Notes on the calculation of the virtual distal metatarsal 2 VDM2 
% (VDM2 is on the line between DM1 and DM5, perpendicular to DM2)
% RVDM2 = DM1 + cos() *||DM2-DM1|| *[baselineRDM] 
% RVDM2 = DM1 + dot(DM2-DM1,DM5-DM1)/(||DM2-DM1||*||DM5-DM1||) * ||DM2-DM1|| *[baselineRDM]
% RVDM2 = DM1 + dot(DM2-DM1,DM5-DM1) *[baselineRDM]
RVDM2  = (Markers.right_DM1_marker + dot((Markers.right_DM2_marker - Markers.right_DM1_marker),baselineRDM)); % right virtual DM2 marker
rmffjc = (Markers.right_PM1_marker + RVPM2)/2; % right medial forefoot joint centre
rlffjc = (Markers.right_PM5_marker + RVPM2)/2; % right lateral forefoot joint centre
rhlxjc = (Markers.right_DM1_marker + RVDM2)/2; % right hallux joint centre
rdigjc = (Markers.right_DM5_marker + RVDM2)/2; % right digiti joint centre

Segments.rHindfoot    = createCS(rhfjc, rmfjc-rhfjc,   Markers.right_CUB_marker - Markers.right_NAV_marker,'zxy');
Segments.rMidfoot     = createCS(rmfjc, rmfdc-rmfjc,   Markers.right_PM5_marker - Markers.right_PM1_marker,'zxy');
Segments.rMedforefoot = createCS(rmffjc,rhlxjc-rmffjc, Markers.right_DM5_marker - Markers.right_DM1_marker,'zxy');
Segments.rLatforefoot = createCS(rlffjc,rdigjc-rlffjc, Markers.right_DM5_marker - Markers.right_DM1_marker,'zxy');
Segments.rHallux      = createCS(rhlxjc,Markers.right_distal_hallux_marker -rhlxjc, Markers.right_DM5_marker- Markers.right_DM1_marker,'zxy');

% left foot segments RDub
lhfjc = (Markers.left_calcaneus_marker+ Markers.left_calcaneus1_marker)/2; % right hindfoot joint centre
lmfjc = (Markers.left_NAV_marker + Markers.left_CUB_marker)/2;  % right midfoot joint centre
lmfdc = (Markers.left_PM1_marker + Markers.left_PM5_marker)/2;  % right midfoot distal centre
% virtual/additional points
LVPM2   = ( 2 * Markers.left_PM1_marker + Markers.left_PM5_marker )/3; % right virtual PM2 marker
baselineLDM = ( Markers.left_DM5_marker - Markers.left_DM1_marker );   % right distal forefoot baseline
baselineLDM = normLength(baselineLDM);
LVDM2  = (Markers.left_DM1_marker + dot((Markers.left_DM2_marker - Markers.left_DM1_marker),baselineLDM)); % right virtual DM2 marker
lmffjc = (Markers.left_PM1_marker + LVPM2)/2; % right medial forefoot joint centre
llffjc = (Markers.left_PM5_marker + LVPM2)/2; % right lateral forefoot joint centre
lhlxjc = (Markers.left_DM1_marker + LVDM2)/2; % right hallux joint centre
ldigjc = (Markers.left_DM5_marker + LVDM2)/2; % right digiti joint centre

Segments.lHindfoot    = createCS(lhfjc, lmfjc -lhfjc,  Markers.left_CUB_marker - Markers.left_NAV_marker,'zxy');
Segments.lMidfoot     = createCS(lmfjc, lmfdc -lmfjc,  Markers.left_PM5_marker - Markers.left_PM1_marker,'zxy');
Segments.lMedforefoot = createCS(lmffjc,lhlxjc-lmffjc, Markers.left_DM5_marker - Markers.left_DM1_marker,'zxy');
Segments.lLatforefoot = createCS(llffjc,ldigjc-llffjc, Markers.left_DM5_marker - Markers.left_DM1_marker,'zxy');
Segments.lHallux      = createCS(lhlxjc,Markers.left_distal_hallux_marker -lhlxjc, Markers.left_DM5_marker - Markers.left_DM1_marker,'zxy');
end
