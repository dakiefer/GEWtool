function progressbar(varargin) % message, step, total, binWidth
% PROGRESSBAR - Display a progress bar in the command window. 
%
% usage: 
% progressbar(message);  % call before the loop with a string to initialize. 
% progressbar(step, total); % current step number and total number of steps
% progressbar(step, total, binWidth); % additionally specify the binwidht in %
%
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    persistent msg;      % message to print: [] by default
    persistent binInd;   % current bin we are at
    if nargin == 1 % reset progressbar
        msg = varargin{1};
        binInd = 0;
        return;
    elseif nargin == 2
        step = varargin{1};
        total = varargin{2}; 
        binWidth = 2; % bins are of width 2 percent per default
    elseif nargin == 3
        step = varargin{1};
        total = varargin{2}; 
        binWidth = varargin{3};
    else
        error('wrong call.')
    end
    nBins = 100/binWidth;
    per = round(step/total*100); 
    if per >= binInd
        nBars = round(per/100*nBins);
        strBars = repmat('❚',1,nBars); % append a bar
        strEmpty = repmat(' ', 1, nBins - nBars);
        lenStr = length(msg) + nBins + 4 + 4;
        strDel = repmat('\b', 1, lenStr); % printing \b moves the cursor back one position
        if binInd ~= 0
            fprintf(strDel);  % if bar has been printed already, go back to beginning.
        end
        fprintf([msg ' [' strBars strEmpty '] ' '%2.0d%%\n'], per); 
        pause(0.01); % updated terminal output
        binInd = binInd + binWidth;
    end
    if step == total  % finished loop
        binInd = 0;   % reset for next use 
    end
end
