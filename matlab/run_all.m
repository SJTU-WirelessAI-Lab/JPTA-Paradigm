function run_all(which)
%RUN_ALL  Reproduce the JPTA 2022 results, one entry point.
%
%   run_all            % self-test + Fig. 4
%   run_all('selftest')
%   run_all('fig4')
%
%   Console output is ASCII on purpose: MATLAB -batch on a Chinese Windows
%   console prints UTF-8 into a GBK code page and the text is lost.

if nargin < 1 || isempty(which)
    which = 'all';
end

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'core'), fullfile(here,'alg'), fullfile(here,'figs'));

t00 = tic;
switch lower(which)
    case 'selftest'
        selftest();
    case 'fig4'
        cfg = cfg_default('N',64,'kappa',64,'maxiter',10);
        fig4_arraygain(cfg);
    case 'all'
        selftest();
        cfg = cfg_default('N',64,'kappa',64,'maxiter',10);
        fig4_arraygain(cfg);
    otherwise
        error('run_all: unknown target "%s"', which);
end
fprintf('\n[run_all] total %.1f s\n', toc(t00));
end
