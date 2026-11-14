function A = array_resp(cfg, theta)
%ARRAY_RESP 式(2) 的阵列响应 a_k(theta)，M x K。
%
%   [a_k(theta)]_m = exp( j*pi*(m-1)*sin(theta)*f_k/f0 ),  m = 1..M
%
%   注意指数从 (m-1) 起，且带 f_k/f0 因子（保留波束斜视）。
%   写成 m*pi*sin(theta) 会让整幅热图斜一格，并与算法 2 的
%   phi_m = pi*(m-1)*sin(theta0) 不自洽。
%
%   theta 为标量时所有子载波用同一角度；为 K x 1 向量时第 k 列用 theta(k)
%   （行为 1 的期望波束就是这么定义的）。

m = (0:cfg.M-1).';
fr = cfg.f(:).' / cfg.f0;          % 1 x K

if isscalar(theta)
    A = exp(1j*pi*m*(sin(theta)*fr));
else
    A = exp(1j*pi*m*(sin(theta(:)).'.*fr));
end
end
