function s=lmitbx_options(options)
% function s=lmitbx_options(options)
%
% displays and modifies the options parameter
% for the LMI Control Toolbox (see "help mincx" for details)
% with no argument, displays the present content
%
% options=[accuracy iterations radius progress trace]
%
% accuracy: relative optimization error
% iterations: maximal number of
% radius: feasibility radius
% progress: slow progress iterations number
% trace: turn-off the trace
% 
% Written by ameg@mit.edu, last modified October 18, 1997

% safety checks
global ABST
if ~isfield(ABST,'options'), 
   error('options field not provided')
end

if nargin>0,
   ABST.options=options;
end

options=ABST.options;
if options(1)==0, options(1)=0.01; end
if options(2)==0, options(2)=100; end
if options(3)==0, options(3)=1e9; end
if options(4)==0, options(4)=10; end

%     12345678901234567890123456789012345678901234567890
s=[  'LMI Control Toolbox "mincx" options:              '];
s=[s;'**************************************************'];
s=[s;'Relative accuracy:            ' sft(options(1),20)];
s=[s;'Maximal # of iterations:      ' sft(options(2),20)];
s=[s;'Feasibility radius:           ' sft(options(3),20)];
s=[s;'Slow progress #:              ' sft(options(4),20)];
s=[s;'Turn-off trace:               ' sft(options(5),20)];
