function [b, bbar, absalpha] = gen_beam(cfg, behavior)
%GEN_BEAM  生成期望波束 {b_k}（第三节），及其归一化版本与 |alpha_k|。
%
%   behavior 1: b_k = sqrt(Psum/(M*K)) * a_k( theta0 + k*dtheta/K )
%   behavior 2: b_k = sqrt(Psum/(M*K)) * a_k(theta1)   , k <  0
%               b_k = sqrt(Psum/(M*K)) * a_k(theta2)   , k >= 0
%
%   输出：
%     b        M x K  期望波束（含功率分配）
%     bbar     M x K  归一化期望波束，每列单位范数
%     absalpha K x 1  |alpha_k| = ||b_k|| = sqrt(Psum/K)

if nargin < 2 || isempty(behavior)
    behavior = cfg.behavior;
end

switch behavior
    case 1
        th = cfg.theta0 + cfg.k*cfg.dtheta/cfg.K;
        A  = array_resp(cfg, th(:));
    case 2
        th = cfg.theta2*ones(cfg.K,1);
        th(cfg.k < 0) = cfg.theta1;
        A  = array_resp(cfg, th(:));
    otherwise
        error('gen_beam: 未知的 behavior=%d', behavior);
end

b        = sqrt(cfg.Psum/(cfg.M*cfg.K)) * A;
bbar     = A/sqrt(cfg.M);                       % 每列范数为 1
absalpha = sqrt(cfg.Psum/cfg.K)*ones(cfg.K,1);  % 式(5) 的最优解
end
