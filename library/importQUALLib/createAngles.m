function out = createAngles(csA,csB,order,Label)
%% create the kardan sequence for the angle between two local coordinate systems (cs) of size[4 3 ns]
% csA   : [4  3 nsamp] array of matrices, forming a time series of a local coordinate sysem
% csB   : [4  3 nsamp] array of matrices, forming a time series of a local coordinate sysem
% order : Kardan sequence
% out   : [3 nsamp] array of Kardan angle
%
% (c) 2020 by Predimo GmbH
% 221004 (MdL) FIXED : handle data with nan values

matvecA = csA(1:3,:,:); % omit the position from the coordinate system size[4 3 ns]
matvecB = csB(1:3,:,:); %
matvecB = permute(matvecB,[2 1 3]); % transpose the matrixB

% do multiplication of the arrays of matrices
ABvec = multiplyMatrixArrays(matvecA,matvecB);

% compute the angles, depending on the Kardan order
if     strcmp(order,'xyz')
    thy = squeeze( asin( ABvec(1,3,:)));
    thz = atan2(-squeeze(ABvec(1,2,:)), squeeze(ABvec(1,1,:)));
    thx = atan2(-squeeze(ABvec(2,3,:)), squeeze(ABvec(3,3,:)));
    out = [thx thy thz]';
elseif strcmp(order,'zyx')
    thy = squeeze( asin(-ABvec(3,1,:)));
    thz = atan2( squeeze(ABvec(2,1,:)), squeeze(ABvec(1,1,:))); 
    thx = atan2( squeeze(ABvec(3,2,:)), squeeze(ABvec(3,3,:)));
    out = [thz thy thx]';
elseif strcmp(order,'zxy')
    thx = squeeze( asin( ABvec(3,2,:)));
    thy = atan2(-squeeze(ABvec(3,1,:)), squeeze(ABvec(3,3,:)));
    thz = atan2(-squeeze(ABvec(1,2,:)), squeeze(ABvec(2,2,:)));
    out = [thz thx thy]';
elseif strcmp(order,'yxz')
    thx = squeeze( asin(-ABvec(2,3,:)));
    thy = atan2( squeeze(ABvec(1,3,:)), squeeze(ABvec(3,3,:)));
    thz = atan2( squeeze(ABvec(2,1,:)), squeeze(ABvec(2,2,:)));
    out = [thy thx thz]';
elseif strcmp(order,'yzx')
    thz = squeeze( asin( ABvec(2,1,:)));
    thx = atan2(-squeeze(ABvec(2,3,:)), squeeze(ABvec(2,2,:)));
    thy = atan2(-squeeze(ABvec(3,1,:)), squeeze(ABvec(1,1,:)));
    out = [thy thz thx];
elseif strcmp(order,'xzy')
    thz = squeeze( asin(-ABvec(1,2,:)));
    thy = atan2(squeeze( ABvec(1,3,:)), squeeze(ABvec(1,1,:)));
    thx = atan2(squeeze( ABvec(3,2,:)), squeeze(ABvec(2,2,:)));
    out = [thx thz thy]';
end

%% Unwrap: detect and repair flips crossing the -180/180 degrees    % MdL 20200607
out      = unwrap(out,[],2);

%% Handle Gimbal locks, which occur in extended position of joints. % MdL 20200602
% -> first and third axis are aligned, and about zero with combined pi-rad flippings
% detect the flips. Take factor 0.99 because the jump can be slightly be larger and smaller than pi
flips    = [[0; 0; 0] round(0.99*diff(out,1,2)/pi)];
flips(isnan(flips)) = 0;       % the data may contain nan values, and at these locations there are no flips
cumflips = pi*cumsum(flips,2); % add or subtract multiples of pi after each flip
out2     = out - cumflips;     % by subtraction, all the flips are removed "unwrap with pi rather than 2*pi"

% The flips in the 1st and 3rd axis cause reversal in the angular velocity in the 2nd axis.
g1=order=='x'; g2=find(order=='y'); g3=order=='z'; % determine the gimbal order
Reversal  = abs(flips(g1,:) .* flips(g3,:)); % detect reversals of both 1st and 3rd axis
Invert    = ones(1,length(Reversal));      % create array of 1 and -1 to invert the angle
CumInvert = cumsum([0 Reversal(1:end-1)]); % increases by 1 at each flip
Invert(mod(CumInvert,2)==1) = -1;          % make -1 at the odd intervals
Idx       = find(Reversal);             % the indices of the reversals
Val       = out(2,Idx);                 % the angular value at the reversal in the 2nd axis
Add       = zeros(1,length(Reversal));  % create Add, the initial angula value of the 2nd axis in the odd periods
for i=2:2:length(Idx)
    for j=Idx(i-1)+1 : Idx(i)
        Add(j) = Val(i-1);
    end
end
% At the odd periods, the 2nd axis is multiplied by -1 and twice the initial value is added
out2(g2,:) = out2(g2,:) .* Invert + 2 * Add;

% Now, a gimbal lock can occur in the g2 component (example DHB_Vis_0094 of sept 2019)   % bug 20200615-MdL
flips    = [[0; 0; 0] round(0.99*diff(out2,1,2)/pi)];
flips(isnan(flips)) = 0;       % the data may contain nan values, and at these locations there are no flips
cumflips = pi*cumsum(flips,2); % add or subtract multiples of pi after each flip
out2     = out2 - cumflips;     % by subtraction, all the flips are removed "unwrap with pi rather than 2*pi"

%% plotting if desired
if 0==1 && (strcmp('rKnee',Label) || strcmp('lKnee',Label))
    figure; title('create\_angles : Check flipping of the circularity crossings')
    subplot(1,3,1); hold on; plot(out(1,:));plot(out2(1,:)); yline(-pi); yline(pi); title([order(1) Label]);
    subplot(1,3,2); hold on; plot(out(2,:));plot(out2(2,:)); yline(-pi); yline(pi); title([order(2) Label]);
    subplot(1,3,3); hold on; plot(out(3,:));plot(out2(3,:)); yline(-pi); yline(pi); title([order(3) Label]);
end

%convert to degrees
out = out2*180/pi;
end
