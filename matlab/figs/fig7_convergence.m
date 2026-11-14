function fig7_convergence(cfg0, nreal, recompute)
%FIG7_CONVERGENCE  图 7：算法 1 的收敛速度（分位带）。
%
%   fig7_convergence()                  % 默认 100 组，有缓存就只重绘
%   fig7_convergence(cfg0, 100, true)   % 强制重算
%
%   定义（第七节）：
%     Fobj(i) = 取 max-iter = i 时算法 1 达到的 Fobj
%     纵轴    = Fobj(i)/Fobj(30)
%     统计对象 = "many realizations" 的 (kappa, N, theta0, Dtheta)（行为1）
%                与 (kappa, N, theta1, theta2)（行为2）
%   画 mean / 10-percentile / 90-percentile，并在 10–90 之间填充"Feasible region"。
%
%   论文未给随机参数的具体分布，本实现取（与 REPORT.zh.md 第 3.4 节一致）：
%     N ∈ {1,2,4,8,16,32,64}（必须整除 M）
%     kappa ∈ [1, M]
%     行为1: theta0 ∈ [-pi/3, pi/3]，Dtheta ∈ [0, pi/2]
%     行为2: theta1, theta2 ∈ [-pi/3, pi/3]
%   随机种子固定为 1，保证可复算。

if nargin < 1, cfg0 = []; end
if nargin < 2 || isempty(nreal), nreal = 100; end
if nargin < 3 || isempty(recompute), recompute = false; end
if isempty(cfg0), cfg0 = cfg_default(); end

outdir  = fullfile(fileparts(mfilename('fullpath')),'..','out');
matfile = fullfile(outdir, sprintf('fig7_convergence_n%d.mat', nreal));
ITMAX   = 30;                       % 归一化基准 Fobj(30)

if exist(matfile,'file') && ~recompute
    S = load(matfile);
    iters = S.iters; Q = S.Q;
    fprintf('[fig7] loaded cached data (%s)\n', matfile);
else
    nf = [1 2 4 8 16 32 64];        % N 的候选（须整除 M）
    rng(1);
    iters = (1:ITMAX).';
    Q = zeros(ITMAX,3,2);           % [i x {p10, mean, p90} x behavior]
    Ns = zeros(nreal,2);            % 记录用过的 N，便于复核
    t0 = tic;

    for beh = 1:2
        R = zeros(ITMAX, nreal);    % 每个 realization 一列
        for r = 1:nreal
            Nn  = nf(randi(numel(nf)));
            kap = 1 + (cfg0.M-1)*rand;
            args = {'N',Nn, 'kappa',kap, 'maxiter',ITMAX};
            if beh == 1
                th0 = -pi/3 + (2*pi/3)*rand;
                dth = (pi/2)*rand;
                args = [args, {'theta0',th0, 'dtheta',dth}]; %#ok<AGROW>
            else
                args = [args, {'theta1', -pi/3 + (2*pi/3)*rand, ...
                               'theta2', -pi/3 + (2*pi/3)*rand}]; %#ok<AGROW>
            end
            cfg = cfg_default(args{:});
            [~, bbar] = gen_beam(cfg, beh);
            [~, ~, ~, hist] = alg1_jpta(cfg, bbar, 'ls');
            R(:,r) = hist/hist(end);
            Ns(r,beh) = Nn;

            if mod(r,10)==0
                fprintf('  beh%d %3d/%3d  (%.1fs)\n', beh, r, nreal, toc(t0));
            end
        end
        for i = 1:ITMAX
            Q(i,1,beh) = pctile(R(i,:), 10);
            Q(i,2,beh) = mean(R(i,:));
            Q(i,3,beh) = pctile(R(i,:), 90);
        end
    end
    save(matfile, 'iters','Q','Ns','-v7');
end

% ------------------------------------------------------------------ 绘图
xr = {[1 ITMAX], [1 15]};
yr = {[0.4 1.0], [0.85 1.0]};
caps = {'(a) Behavior 1','(b) Behavior 2'};

fig = figure('Visible','off','Position',[100 80 620 860],'Color','w');

for beh = 1:2
    ax = subplot(2,1,beh);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    p10 = Q(:,1,beh); mu = Q(:,2,beh); p90 = Q(:,3,beh);

    % Feasible region：10–90 分位之间
    h0 = fill(ax, [iters; flipud(iters)], [p10; flipud(p90)], ...
              [1 1 0], 'EdgeColor','none');
    h1 = plot(ax, iters, p10, 'k:',  'LineWidth',1.8);
    h2 = plot(ax, iters, mu,  'r-',  'LineWidth',2.2);
    h3 = plot(ax, iters, p90, 'k--', 'LineWidth',1.8);

    xlim(ax, xr{beh}); ylim(ax, yr{beh});
    set(ax,'FontSize',11,'TickDir','out','Layer','top');
    xlabel(ax,'max-iter (i)','FontSize',12);
    ylabel(ax,'F_{obj}(i) / F_{obj}(30)','FontSize',12);

    lg = legend(ax, [h0 h1 h2 h3], ...
                {'Feasible region','10-percentile','Mean','90-percentile'}, ...
                'Location','east','FontSize',10);
    lg.Box = 'on';

    text(ax, 0.5, -0.19, caps{beh}, 'Units','normalized', ...
         'HorizontalAlignment','center','FontSize',12);
end

print(fig, fullfile(outdir,'fig7_convergence.png'), '-dpng','-r200');
close(fig);
fprintf('[fig7] figure written\n');
end

% ----------------------------------------------------------------------
function q = pctile(x, p)
%PCTILE  线性插值分位数（不依赖 Statistics Toolbox）。
y   = sort(x(:));
n   = numel(y);
pos = p/100*(n-1) + 1;
lo  = floor(pos); hi = ceil(pos);
q   = y(lo) + (pos-lo)*(y(hi)-y(lo));
end
