function out = fig4_arraygain(cfg)
%FIG4_ARRAYGAIN  Fig. 4: array gain |a_k(theta)' * T*P*d_k|^2 for the ideal and
%   the JPTA beamformer, beam behaviors 1 and 2.
%
%   Panels: (a) Ideal, Behavior 1   (b) JPTA, Behavior 1
%           (c) Ideal, Behavior 2   (d) JPTA, Behavior 2
%
%   Paper setting: N = M = 64, kappa = M, max-iter = 10 for Algorithm 1.
%   Layout copied from the original: axis labels, x ticks 0.95/1/1.05 in
%   units of 1e11 Hz, one "Array Gain" colorbar per panel, 0..M linear
%   (the peak array gain is M = 64 -- the scale is NOT in dB).
%   Added on top of the original: the achieved Fobj printed inside (b) and (d).
%
%   Usage: cfg = cfg_default('N',64,'kappa',64); out = fig4_arraygain(cfg);

if nargin < 1 || isempty(cfg)
    cfg = cfg_default();
end

here = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(here,'core'), fullfile(here,'alg'));
outdir = fullfile(here,'out');
if ~exist(outdir,'dir'), mkdir(outdir); end

theta = (-90:0.25:90).';
fprintf('Fig.4: array gain on %d angles x %d subcarriers ...\n', ...
        numel(theta), cfg.K);

G = cell(1,4); tau = cell(1,2); phi = cell(1,2); F = zeros(1,2);

for b = 1:2
    [~, bbar] = gen_beam(cfg, b);
    G{2*b-1} = array_gain(cfg, bbar, theta);            % ideal
    t0 = tic;
    [tau{b}, phi{b}] = alg1_jpta(cfg, bbar, 'ls');      % JPTA, eq.(8a)
    X = jpta_forward(cfg, tau{b}, phi{b});
    G{2*b} = array_gain(cfg, X, theta);
    F(b)   = fobj(cfg, bbar, X);
    fprintf('  behavior %d : Fobj = %.4f  (%.1f s)\n', b, F(b), toc(t0));
end

% ---- plot -------------------------------------------------------------
caps  = {'(a) Ideal, Behavior 1', '(b) JPTA, Behavior 1', ...
         '(c) Ideal, Behavior 2', '(d) JPTA, Behavior 2'};
fnote = {'', sprintf('F_{obj} = %.3f', F(1)), ...
         '', sprintf('F_{obj} = %.3f', F(2))};

fig = figure('Visible','off','Position',[100 80 960 900],'Color','w');

xf = cfg.f/1e11;                       % paper: f in units of 1e11 Hz
for i = 1:4
    ax = subplot(2,2,i);
    imagesc(ax, xf, theta, G{i});
    set(ax,'YDir','normal','FontSize',10,'TickDir','out','Box','on');
    axis(ax,'xy');
    colormap(ax, parula);
    clim(ax, [0 cfg.M]);
    xlim(ax, [xf(1) xf(end)]);
    ylim(ax, [-75 75]);
    xticks(ax, [0.95 1 1.05]);
    yticks(ax, [-50 0 50]);
    % 略微压低坐标区，给下方的面板标题腾位置
    p = ax.Position; ax.Position = [p(1) p(2)+0.03 p(3) p(4)*0.93];
    cb = colorbar(ax);
    cb.FontSize = 9; cb.Title.String = 'Array Gain'; cb.Title.FontSize = 9;
    xlabel(ax,'Frequency (Hz) \times10^{11}','FontSize',10);
    ylabel(ax,'Angle of Departure \theta (deg)','FontSize',10);
    % panel caption below the axes, as in the original
    text(ax, 0.5, -0.30, caps{i}, 'Units','normalized', ...
         'HorizontalAlignment','center','FontSize',13,'Interpreter','none');
    if ~isempty(fnote{i})
        text(ax, 0.03, 0.94, fnote{i}, 'Units','normalized', ...
             'Color','w','FontWeight','bold','FontSize',10, ...
             'BackgroundColor',[0 0 0 0.35],'Margin',3);
    end
end

png = fullfile(outdir,'fig4_arraygain.png');
exportgraphics(fig, png, 'Resolution', 150);
close(fig);

mat = fullfile(outdir,'fig4_arraygain.mat');
theta_deg = theta; f_Hz = cfg.f;
G4 = cellfun(@single, G, 'UniformOutput', false);   % single 存盘，体积减半
tau4 = tau; phi4 = phi; F4 = F; cfg4 = cfg;
save(mat,'theta_deg','f_Hz','G4','tau4','phi4','F4','cfg4','-v7.3');

fprintf('  saved %s\n  saved %s\n', png, mat);

if nargout > 0
    out = struct('png',png,'mat',mat,'Fobj',F,'G',{G},'theta',theta);
end
end
