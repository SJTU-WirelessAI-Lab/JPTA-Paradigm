function fig9_hbf_arraygain(cfg, recompute)
%FIG9_HBF_ARRAYGAIN  图 9：HBF-FC 与 HBF-PC 的阵列增益热图。
%
%   fig9_hbf_arraygain()             % 有缓存就只重绘
%   fig9_hbf_arraygain(cfg, true)    % 强制重算
%
%   四个面板（照原文）：
%     (a) HBF-FC, N_RF=22, 行为 1   色标 0-70
%     (b) HBF-FC, N_RF= 2, 行为 2   色标 0-60
%     (c) HBF-PC, N_RF=32, 行为 1   色标 0-70
%     (d) HBF-PC, N_RF=32, 行为 2   色标 0-60
%
%   增益定义与图 4 相同：|a_k(theta)' * x_k|^2，x_k = F_RF*f_BB,k 单位化。
%   对比图 4 可见 HBF 的波束指向随频率的跳变更"锐"，这是论文的论点。

if nargin < 1 || isempty(cfg), cfg = cfg_default(); end
if nargin < 2 || isempty(recompute), recompute = false; end

outdir  = fullfile(fileparts(mfilename('fullpath')),'..','out');
matfile = fullfile(outdir,'fig9_hbf_arraygain.mat');

spec = struct( ...
    'name', {'(a) HBF-FC, N_{RF} = 22, Behavior 1', ...
             '(b) HBF-FC, N_{RF} = 2, Behavior 2', ...
             '(c) HBF-PC, N_{RF} = 32, Behavior 1', ...
             '(d) HBF-PC, N_{RF} = 32, Behavior 2'}, ...
    'kind', {'fc','fc','pc','pc'}, ...
    'nrf',  {22, 2, 32, 32}, ...
    'cmax', {70, 60, 70, 60}, ...
    'ylim', {[-75 75], [-75 75], [-75 75], [-90 90]}, ...
    'yticks', {[-50 0 50], [-50 0 50], [-50 0 50], [-90 -45 0 45 90]});
behlist = [1 2 1 2];

theta = -90:0.5:90;

if exist(matfile,'file') && ~recompute
    S = load(matfile);
    G9 = S.G9; theta = S.theta;
    fprintf('[fig9] loaded cached data (%s)\n', matfile);
else
    G9 = cell(1,4);
    t0 = tic;
    for i = 1:4
        b = gen_beam(cfg, behlist(i));
        [~, bbar] = gen_beam(cfg, behlist(i));
        if strcmp(spec(i).kind,'fc')
            [F, Bb] = hbf_fc(b, spec(i).nrf, [], [], bbar);
        else
            [F, Bb] = hbf_pc(b, spec(i).nrf, [], [], bbar);
        end
        X = F*Bb; X = X./vecnorm(X);
        G9{i} = array_gain(cfg, X, theta);
        fprintf('  %s : done (%.1fs)\n', spec(i).name, toc(t0));
    end
    save(matfile, 'G9','theta','-v7');
end

% ------------------------------------------------------------------ 绘图
fig = figure('Visible','off','Position',[100 60 980 900],'Color','w');
xf  = cfg.f/1e11;

for i = 1:4
    ax = subplot(2,2,i);
    imagesc(ax, xf, theta, G9{i});
    set(ax,'YDir','normal','FontSize',10,'TickDir','out','Box','on');
    axis(ax,'xy');
    colormap(ax, parula);
    clim(ax, [0 spec(i).cmax]);
    xlim(ax, [xf(1)-5e-4, xf(end)+5e-4]);   % 留一点余量，否则 1.05 的刻度不显示
    ylim(ax, spec(i).ylim);
    xticks(ax, [0.95 1 1.05]);
    yticks(ax, spec(i).yticks);
    cb = colorbar(ax);
    cb.FontSize = 9; cb.Title.String = 'Array Gain'; cb.Title.FontSize = 9;
    xlabel(ax,'Frequency (Hz) \times10^{11}','FontSize',10);
    ylabel(ax,'Angle of Departure \theta (deg)','FontSize',10);
    text(ax, 0.5, -0.27, spec(i).name, 'Units','normalized', ...
         'HorizontalAlignment','center','FontSize',12,'Interpreter','tex');
end

print(fig, fullfile(outdir,'fig9_hbf_arraygain.png'), '-dpng','-r200');
close(fig);
fprintf('[fig9] figure written\n');
end
