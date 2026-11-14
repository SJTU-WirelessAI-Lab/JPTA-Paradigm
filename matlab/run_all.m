function run_all(which)
%RUN_ALL  Reproduce the JPTA 2022 results, one entry point.
%
%   run_all            % self-test + Fig. 4 + Fig. 5 + Fig. 6
%   run_all('selftest')
%   run_all('fig4')
%   run_all('fig5')    % Fobj vs N   (~2-4 min)
%   run_all('fig6')    % Fobj vs kappa (~2-4 min)
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
    case 'fig5'
        fig5_fobj_vs_N(cfg_default());
    case 'fig6'
        fig6_fobj_vs_kappa(cfg_default());
    case 'all'
        selftest();
        cfg = cfg_default('N',64,'kappa',64,'maxiter',10);
        fig4_arraygain(cfg);
        fig5_fobj_vs_N(cfg_default());
        fig6_fobj_vs_kappa(cfg_default());
    otherwise
        error('run_all: unknown target "%s"', which);
end
fprintf('\n[run_all] total %.1f s\n', toc(t00));
end
