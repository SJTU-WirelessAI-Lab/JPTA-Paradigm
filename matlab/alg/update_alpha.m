function ang_alpha = update_alpha(cfg, bbar, tau, phi)
%UPDATE_ALPHA  式(11)：给定 tau、phi 后的 digital 相位。
%
%   ang_alpha_k = angle( sum_n sum_{m in M_n} [bbar_k]_m
%                        * e^{-j*phi_m} * e^{j*2*pi*f_k*tau_n} )

E  = exp(1j*2*pi*(cfg.f(:).').*tau(cfg.gidx));      % M x K，按天线所属 TTD
ang_alpha = angle(sum((bbar.*exp(-1j*phi(:))).*E, 1)).';
end
