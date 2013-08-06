function disp_str(varargin)

warning('on','all')

if nargin < 1
    disp_str(1)
end

idx=varargin{1};
n_parameter=length(varargin)-1;

for i1=1:n_parameter
    parameter{i1}=varargin{i1+1}; %#ok<*AGROW>
end

switch idx
    case 1
        error('input must be specified')
    case 2
        error('input argument too long')
    case 3
        error('output must be specified')
    case 4
        error([parameter{1},' first arguments must be defined'])
        % parameter{1} = 'Two' or 'Three'
    case 5
        error(['Improper ',parameter{1},' argument size'])
        % parameter{1} = 'second' or 'third'
    case 6
        error('Incompatible horizontal size')
    case 7
        error(['input of ''',parameter{1},''' is not greater than ',...
            parameter{2}])
        % parameter{1} = 'skew'
        % parameter{2} = 'one'
        
        
    case 11
        warning('The "abst" environment not defined')
    case 12
        error('"abst" environment not initialized')
    case 13
        error('"abst" solution not available')
    case 14
        error(['Operation "',parameter{1},'" not supported'])
        % parameter{1} = 'and' or 'ctranspose'
    case 15
        error(['This is not an "',parameter{1},'" environment'])
        % parameter{1} = 'iqc'
        
        
        
    case 31
        error('Bad class name')
    case 32
        disp('There is no log entry for the specified class')
    case 33
        error(['Bad index! Log has only #',parameter{1},' entries.'])
        % parameter{1} = num2str(ABST.nlog)
    case 34
        error('argument of "trace" must be square')
    case 35
        error('The log entry referenced here does not exist')
    case 36
        error(['***',parameter{1},': conversion ',parameter{2},...
            '->abst  not allowed'])
        % parameter{1} = ABST.name
        % parameter{2} = intype
    case 37
        error('Improper use of "abst" constructor')
    case 38
        error('argument dimensions are not compatible')
    case 39
        error([parameter{1},' ',parameter{2},' ',parameter{3},...
            ' not allowed'])
        % parameter{1} = ABST.cls{ca}
        % parameter{2} = '&'
        % parameter{3} = ABST.cls{cb}
    case 40
        error([parameter{1},parameter{2},'  not allowed'])
        % parameter{1} = ABST.cls{ca}
        % parameter{2} = ''''
    case 41
        error('only a signal or an input could be differentiated')
    case 42
        error(['[',parameter{1},parameter{2},parameter{3},...
            '] not allowed'])
        % parameter{1} = ABST.cls{ca}
        % parameter{2} = ' '
        % parameter{3} = ABST.cls{cb}
    case 43
        error([parameter{1},parameter{2},parameter{3},' not allowed'])
        % parameter{1} = ABST.cls{ca}
        % parameter{2} = [s]=
        % parameter{3} = ABST.cls{cb}
    case 44
        warning('the "LMI tool" is not available for the user')
    case 45
        warning(['the "',parameter{1},...
            '" is not available for the user(''lmilab'' or ''yalmip'')'])
        % parameter{1} = lmitool
    case 46
        warning(['Parameter setting error! The parameter matrix',...
            ' must be a 1x5 vector; e.g. [0 200 0 0 0].'])
    case 47
        warning(['Parameter setting error! The parameter has ',...
            'the form of ''parameter name'',value,',...
            '''parameter name'',value,....'])
    case 48
        disp(['The ',parameter{1},' default value is used!'])
        % parameter{1} = lmilab or yalmip
    case 49
        error(['You must run an optimizer (e.g. "',...
            parameter{1},'") first'])
        % parameter{1} = iqc_gain_tbx
    case 50
        disp([parameter{1},'(exporting the LMI script to ',parameter{1},...
            '_exe.m) ...'])
        % parameter{1} = iqc_analysis_l2gain_lmilab
    case 51
        disp([parameter{1},' ...'])
        % parameter{1} = iqc_analysis_l2gain_lmilab
    case 52
        disp('  defining the original variables ...')
    case 53
        disp('  defining the non-KYP LMIs ...')
    case 54
        disp('  defining the KYP LMIs ...')
    case 55
        error('Unknown error...')
    case 56
        disp(['  Solving with ',parameter{1},' decision variables ...'])
        % parameter{1} = num2str(ndec)
    case 57
        disp(['Minimize ',parameter{1},' LMI is ',parameter{2},'!'])
        % parameter{1} = L2-gain
        % parameter{2} = Feasible or Infeasible
    case 58
        warning(['Numerical problems, the smallest eigenvalue is ',...
            parameter{1}])
        % parameter{1} = num2str(min(alleig))
    case 59
        warning(['Something else happened: ',parameter{1}]);
        % parameter{1} = yalmiperror(sol.problem)
    case 60
        error('all inputs must be from the "abst" class')
    case 61
        error([parameter{1},' argument must be input or signal'])
        % parameter{1} = First or Second
    case 62
        error('Bad cost function, not a ''var'' or a ''lin''')
    case 63
        error('cost function need to have size 1x1')
    case 64
        disp([parameter{1},': processing the abst log information ... '])
        % parameter{1} = iqc_extract
    case 65
        disp(['  ',parameter{1},' done OK'])
        % parameter{1} = fdlmi_extract
    case 66
        error('Bad index of independent variables!!')
    case 67
         disp(['The problem is ',parameter{1}])
         % parameter{1} = feasible or infeasible
    case 68
         warning(['range auto set [',parameter{1},' ',parameter{2},']!'])
         % parameter{1} = 0
         % parameter{2} = pi
end

warning('off','all')
