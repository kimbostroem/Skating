% reset Matlab to its initial state

% close model, if necessary
if slreportgen.utils.isModelLoaded('Myonardo')
    set_param('Myonardo','FastRestart',0);
    close_system('Myonardo',0);
end

close all; % close open figures
delete(findall(0,'type','figure')); % delete also hidden figures
clear; % clear workspace
clear loadBody; % clear persistent variables in function
restoredefaultpath; % restore default path
rehash toolboxcache % rehash toolbox cache
matlabrc; % reset MATLAB's internal state
clc; % clear command window