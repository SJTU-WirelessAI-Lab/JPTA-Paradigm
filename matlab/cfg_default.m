function cfg = cfg_default(varargin)
%CFG_DEFAULT  JPTA 复现的统一参数表（论文第七节）。
%
%   所有仿真参数集中在这里，图脚本只覆盖 1~2 个字段：
%       cfg = cfg_default('N', 16, 'kappa', 32);
%
%   论文：V. V. Ratnam, J. Mo, A. Alammouri, B. L. Ng, J. Zhang, A. F. Molisch,
%   "Joint Phase-Time Arrays: A Paradigm for Frequency-Dependent Analog
%   Beamforming in 6G", IEEE Access, vol. 10, pp. 73364-73377, 2022.

% ---- 频域（第七节）----------------------------------------------------
cfg.f0   = 100e9;                 % 载频
cfg.W    = 10e9;                  % 系统带宽（定义用，见下）
cfg.K    = 2048;                  % 子载波数
cfg.k    = (-cfg.K/2 : cfg.K/2-1).';                 % -1024 ... 1023
cfg.f    = cfg.f0 + cfg.W*cfg.k/cfg.K;               % 95 ... 104.995 GHz

% ---- 阵列（第三节）----------------------------------------------------
cfg.M    = 64;                    % 天线数，f0 处半波长间距
cfg.N    = 64;                    % TTD 数，N <= M
cfg.kappa= 64;                    % TTD 时延范围 [0, kappa/W]

% ---- 波束行为（第七节）------------------------------------------------
cfg.behavior = 1;
cfg.theta0   = pi/6;              % 行为 1 扫描中心
cfg.dtheta   = pi/4;              % 行为 1 扫描范围
cfg.theta1   = -pi/4;             % 行为 2 低频段角度
cfg.theta2   = pi/6;              % 行为 2 高频段角度

% ---- 算法 -------------------------------------------------------------
cfg.maxiter  = 10;                % 算法 1 迭代数（第七节）
cfg.method   = 'ls';              % 'ls' = 式(8a) 线搜索；'lemma1' = 式(9)(10)
cfg.Nfft     = 32768;             % 式(8a) 的 τ 网格：FFT 补零长度

% ---- 功率与加权 -------------------------------------------------------
cfg.Psum  = 1;
cfg.omega = ones(cfg.K,1);        % ω_k = 1（第七节）

% ---- 用户覆盖 ---------------------------------------------------------
if mod(numel(varargin),2) ~= 0
    error('cfg_default: 覆盖参数必须成对出现 (name, value)');
end
for i = 1:2:numel(varargin)
    cfg.(varargin{i}) = varargin{i+1};
end

% ---- 派生量（必须在覆盖之后重算）--------------------------------------
cfg.fmax  = max(cfg.f);
cfg.fmin  = min(cfg.f);
cfg.tau_max = cfg.kappa/cfg.W;
cfg.P     = build_map(cfg.M, cfg.N);   % 式(1) 的 P：M x N，每行一个 1
cfg.gidx  = cfg.P*(1:cfg.N).';         % 每根天线所属 TTD 编号，M x 1
end
