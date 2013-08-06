function lmitbx_options(options)
% function lmitbx_options(options)
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
% Last modified by cmj on 2013/4/23

% safety checks
global ABST
if ~isfield(ABST,'log'),
    disp_str(12)
end

if nargin==0
    switch ABST.lmitool
        case 'yalmip'
            eval(['setlmioptions(''yalmip'',',ABST.lmiparameter,')'])
        case 'lmilab'
            setlmioptions('lmilab',ABST.lmiparameter)
    end
else
    if options(1)==0, options(1)=0.01; end
    if options(2)==0, options(2)=100; end
    if options(3)==0, options(3)=1e9; end
    if options(4)==0, options(4)=10; end
    setlmioptions('lmilab',options)
end
