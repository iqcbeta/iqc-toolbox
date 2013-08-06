function blkStruct = slblocks
%SLBLOCKS Defines the block library for a specific Toolbox or Blockset.

% Name of the subsystem which will show up in the SIMULINK Blocksets
% and Toolboxes subsystem.
% Example:  blkStruct.Name = 'DSP Blockset';
blkStruct.Name = 'IQC Toolbox';

% The function that will be called when the user double-clicks on
% this icon.
% Example:  blkStruct.OpenFcn = 'dsplib';
blkStruct.OpenFcn = 'iqc_lib';

% The argument to be set as the Mask Display for the subsystem.  You
% may comment this line out if no specific mask is desired.
% Example:  blkStruct.MaskDisplay = 'plot([0:2*pi],sin([0:2*pi]));';
% No display for now.
% blkStruct.MaskDisplay = '';

% Define the library list for the Simulink Library browser.
% Return the name of the library model and the name for it

Browser(1).Library = 'Common';
Browser(1).Name    = 'Common';
Browser(1).IsFlat  = 0;

Browser(2).Library = 'Calculate';
Browser(2).Name    = 'Calculate';
Browser(2).IsFlat  = 0;

Browser(3).Library = 'Simulation';
Browser(3).Name    = 'Simulation';
Browser(3).IsFlat  = 0;

% End of blocks


