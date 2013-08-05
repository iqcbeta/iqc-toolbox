% IQC_beta, version 028, June 01, 2001
%
% Maintained by  Alex     :  ameg@mit.edu 
%                Isaac    :  cykao@mit.edu
% Significant contributions by:
%                Anders   :  rantzer@control.lth.se
%                Fernando :  fdamato@lids.mit.edu
%                Jorge    :  jmg@cds.caltech.edu or jmg@mit.edu (will expire after 2001) 
%                Ulf      :  ulfj@math.kth.se
% Many thanks (for constantly using and reporting bugs) to :
%                Pablo    :  pablo@cds.caltech.edu
% 
% (persons are listed aphabetically)
%
% CURRENT FIXES:
%
% 06.01.01: 1. Change a feature of (fd)lmi_mincx_tbx.m : for feasibility test, 
%              the programs return a positive number as the 'cost' output while feasible, 
%              and a nagative number while infeasible.
%
% 03.09.01: 1. Fix a bug in lmi_extract.m : the structured variables 
%              were handled in a incorrect way
% 
%           2. Add more applications (analysis tools for relay feedback systems, on/off systems,
%              systems with saturation) done by Jorge. 
%
% 12.29.00: 1. Include several works done by Fernando:
%              new version of 
%               (1) iqc_d_slope_odd.m
%               (2) iqc_d_slope.m
%
%              new m-files:
%               (1) getGH_d_slope.m 
%               (2) getGH_d_slope_odd.m
%
%               They build the multipliers G and H from the results of the
%               optimization.
%
%              one example : example_automatica.m
%
% 01.14.00: 1. Fixed an error in iqc_slowtv.m
%
% 08.30.99: 1. Fixed a bug in iqc_ratelimiter.m (reported by Prof. Weerd)  
% 
% 08.19.99: 1. Include some application programs from Jorge.
%
% 07.14.99: 1. Fixed a bug with LMI-environment (square the constant term appearing in 
%              the cost function) discovered by Pablo.
%           2. The sequence of the arguments to iqc_monotonic.m was changed. 
%              Before it was : [v,op,a,k]. Now the 'op' is moved to the last position:
%              [v,a,op,k] !!!  
%           3. Add one more test file (for iqc_ltiunmod.m) : iqc_ltiunmod_test1.m  
% 
% 07.02.99: 1. Fixed two bugs discovered by Pablo.
%           2. Introduce warning messages when "X {+/-/>/<} a" (or "a {+/-/>/<} X")
%              is discovered, where "X" is square and has size large than 1, and 
%              "a" is a scalar. NOTE that, in IQC_beta, "X {+/-/>/<} a" is intepreted as 
%              "X {+/-/>/<} a*ones(size(X))" !! ( NOT "X {+/-/>/<} a*eye(size(X))" )
%
% 06.22.99: In version 017:
%           1. Fixed the problem caused by a bug in "ss/size.m" of MATLAB version 5.3
%           2. Fixed a bug in "lmi_value.m" (discovered by Jorge)
%           3. More test programs added
%           4. All LP based solvers are removed: 
%              Since LP based solvers are still under developement, we don't think now is
%              a good time to publish them. If you are interested in them, please contact
%              ameg@mit.edu or cykao@mit.edu
%
% 04.30.99: Version 016 was prepared by Isaac.  This is the first time Alex 
%           throws this big job to me, which is a good indication that he is 
%           getting bored with this toolbox and something new is probably 
%           coming up. Here is the major changes in this version:
% 
%           1. Re-arrangement of all files: 
%               Because the number of programs is getting bigger and bigger, 
%                I think it makes sense to make things a bit more ordered. 
%                From now on, the directory 'tools' of previous versions is 
%                 splited to 3 directories:
%                1.) absttools : contains programs that relate to ABST environment
%                   (including constructors, all overloading functions, etc.)  
%                2.) iqctools  : contains programs that relate to IQC applications
%                    including the libeary of iqc descriptions, solvers, etc.
%                3.) lmitools  : contains programs that relate to LMI applications
%               
%                The directory structure now is:
% 
%                iqc_beta -- absttools -- @abst
%                         -- iqctools  -- iqclib
%                         -- lmitools
%                         -- tests
% 
%           2. Reversion of LMI environment: 
%                The constraint on forming block of matrix variables is released. 
%                In this new version, (I hope that) matrix variables with any structure 
%               is possible. For example: it is now possible to form:
%  
%                X=[x1,c1,c2,c3,x2; 
%                   c4,c5,x3 c6,c7],  xi: matrix variables , ci: constant matrices.
% 
%                which is not allowed to do in the prevision versions. 
% 
%           3. New environment -- frequency dependent LMI model: 
%                A new environment in which we can define Frequency Dependent LMI's
%                (FDLMI's) is innovated. The initialization of this environment  
%                is "abst_init_fdlmi". This environment inherits most of features
%                from IQC environment but simplified. Following is a simple example
%                to show how it works and why it is useful:
% 
%                Model approximation problem: 
%                   Get a best approximation (in the sense of H-infinity) of 3/(s+3) 
%                   by using 1, 1/(s+1) and 1/(s+2), i.e., find x0, x1, x2,
%                   such that |3/(s+3)-(x0+x1/(s+1)+x2/(s+2)|_inf --> minimized.
% 
%                    abst_init_fdlmi;
%                    G0=ss(-3,3,1,0);
%                    G1=ss(-1,1,1,0);
%                    G2=ss(-2,1,1,0);
%                    x0=symmetric;
%                    x1=symmetric;
%                    x2=symmetric;
%                    y =symmetric;
%                    y>0;
%                    Ga=x0+x1*G1+x2*G2;
%                    H =G0-Ga;
%                    [y,H;H',y]>0;
%                    fdlmi_mincx_tbx(y)
%                 
%                This environment of course also handles the Non-Frequency Dependent
%                LMI's (NFDLMI). However, I will recommend using NFDLMI environment
%                for the problems contain only NFDLMI's. Because all constants in the
%                FDLMI environment are stored and handled as "ss" objects, and we all
%                know that MATLAB Control Toolbox is not very reliable ...  
% 
%            4. LP based solver :
%                Sasha has been talking about this for a long long time and 
%                eventually we found a time to implement this great idea. Three LP
%                based solver are available now:
%                  "iqc_gain_lp.m"    -- solver for IQC problems
%                  "lmi_mincx_lp.m"   -- solver for Non-frequency dependent LMI problems
%                  "fdlmi_mincx_lp.m" -- solver for Frequency dependent LMI problems
%                For how to use them, read the help files associated with those
%                programs.
%                
%                Those solvers are really proto-types and far away from perfect.
%                We are working on elaborating them in every aspect, including
%                the algorithms. Bugs in those programs are expectable. Any time you 
%                find strange things happen, please report to me. (cykao@mit.edu) 
% 
%
% 12.15.98: lmi solver now accepts a wider selection of block structures
%           (previously using something like [1 1;1 x], where x is a variable, would be
%           impossible); iqc_get_mlmi allows access to the "main LMI" of the
%           IQC feasibility LMI problem; iqc_mfig is a new program for drawing
%           simple block diagrams and exporting them in a LaTeX format
% 
% 10.20.98: new function get_pab provides an easy access to
%           the P,A and B matrices from the total IQC feasibility problem
%
% 09.19.98: fixed a bug in lmi_mincx_tbx
%           new version of iqc_gui
%           unmodeled LTI dynamics block
%
% 06.19.98: new versions of iqc_gain_tbx and iqc_extract allow
%           linking ANY signal to ANY input;
%           Zames-Falb analogs for repeated slope-restricted
%           nonlinearities
%
% 06.17.98: fixed a bug in iqc_mincx_tbx and iqc_feas_tbx - 
%           incorrect handling of transposition

