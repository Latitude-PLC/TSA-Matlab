function ccdc_standalone(task,ntasks)
%function myStandalone(n, nprocs)
% function myStandalone(n, nprocs)
% Purpose: this optionally serves as the front m-file to your app
%          so your app remain unchanged.
%       n: size of simple arithmetic sequence [1+2+3+ . . . +n]
%  nprocs: number of processors [ (matlabpool('local', nprocs) ]
% ===> Rules
% a) This program, myStandalone, must be a function WITH no output
%    >> mystandalone 200 4      % command  syntax
%    >> mystandalone(200,4)     % function syntax
% b) The standalone runs in command syntax only
%    scc1% myStandalone 200 4

% Check input arguments
if nargin ~= 2
    disp(['Expecting 2 but received ' num2str(nargin) ...
          ' input parameters, job aborted.'])
    return
end
    
% In command syntax, both n and nprocs are strings
if ischar(task) == 1
    disp(['task = ',task])
    task = str2double(task);   
end

if ischar(ntasks) == 1
    disp(['ntasks = ',ntasks])
    ntasks = str2double(ntasks);   
end

% ccdc application (can be script or function m-file)
main_ChangePar(task,ntasks)

if isdeployed     % in standalone mode . . .
   exit;
else
   close all
end

end   % end of function
