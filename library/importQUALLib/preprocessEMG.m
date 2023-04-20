function ProcessedEMG = preprocessEMG(Rawdata, MessFreq, HochpassFreq, TiefpassFreq, Order)
%% (c) Marc de Lussanet, WWU Muenster
%% Version 1 (20.12.2018) 
%     langsamer Trend entfernen durch high pass filter
%     Rektifizieren, 
%     glaetten mit tiefpassfilter:
%% Version 2 (15.4.2019) 
	
	%% Default values:
	if ~HochpassFreq, HochpassFreq=20; end % (remove artifacts: De Luca et al 2010)
	if ~TiefpassFreq, TiefpassFreq=50; end % 
	if ~Order       , Order = 2;       end
	if mod(Order,2), error('Order must be even!'); end
	
	%% Processing:
	Order = Order/2; % effective order is doubled due to filtfilt
	ProcessedEMG = gapfilth(TiefpassFreq,MessFreq,Order,...
	                     abs( ...
	                         gapfilth(HochpassFreq,MessFreq,Order,Rawdata,'h') ...
	                        ) ...
	                     ,'l');
								
end
