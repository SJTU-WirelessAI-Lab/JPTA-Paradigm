function [tau, info] = solve_tau_ls(cfg, bbar, ang_alpha)
%SOLVE_TAU_LS  Eq.(8a): per-TTD delay by exhaustive line search over the
%   feasible range, for a given set of digital phases.
%
%   tau_n = argmax_tau  sum_{m in M_n} | sum_k omega_k*e^{j*ang_alpha_k}
%                                        * conj([bbar_k]_m)*e^{-j*2*pi*f_k*tau} |
%
%   SEARCH RANGE (important):
%   Eq.(4) searches -kappa/(2W) <= tau_n <= kappa/(2W), i.e. a SYMMETRIC range.
%   Algorithm 1's "[0, kappa/W]" is only the final hardware range, obtained
%   afterwards by the centering step. Searching [0, kappa/W] here is WRONG for
%   beam behavior 1: the delay it needs is NEGATIVE,
%       tau_m ~ -(m-1)*sin(theta0)/(2*f0) ~ -2.1 ns for m = 64,
%   so the peak lies outside [0, kappa/W] and the search locks onto a sidelobe
%   (Fobj drops from ~0.93 to ~0.67). Lemma 1 hits 0.93 because eq.(10a)
%   already uses the symmetric range -- that is the tell.
%
%   IMPLEMENTATION: the candidates are tau = n*K/(W*Nfft), so
%       e^{-j*2*pi*f_k*tau} = e^{-j*2*pi*f0*tau} * e^{-j*2*pi*(j-1)*n/Nfft}
%   where the first factor is a common phase (irrelevant, the objective takes
%   an absolute value) and the second is one FFT bin. One fft therefore yields
%   all M antennas on all Nfft candidate delays; the grid step is
%   K/(W*Nfft) = 6.25 ps and a parabolic interpolation brings it to ~1 ps.
%   The objective is periodic in tau with period K/W, provided f0*K/W is an
%   integer (true here: 100e9*2048/10e9 = 20480), which is what lets the
%   negative half of the range be read off the same FFT array by wrapping.

M = cfg.M; K = cfg.K; Nfft = cfg.Nfft;

c = (cfg.omega.*exp(1j*ang_alpha)).';    % 1 x K
B = conj(bbar).*c;                       % M x K

T = fft(B, Nfft, 2);                     % M x Nfft
step = K/(cfg.W*Nfft);                   % candidate spacing (s)

Nh    = round(cfg.kappa*Nfft/(2*K));     % half-range in grid points
nlist = (-Nh:Nh).';                      % symmetric candidate indices
idx   = mod(nlist, Nfft) + 1;

A   = abs(T(:,idx)).';                   % (2Nh+1) x M
OBJ = A*cfg.P;                           % (2Nh+1) x N, summed per TTD group

[~, imax] = max(OBJ, [], 1);
tau = zeros(cfg.N,1);
for n = 1:cfg.N
    i = imax(n);
    d = 0;
    if i > 1 && i < numel(nlist)         % parabolic refinement
        y   = OBJ(i-1:i+1, n);
        den = y(1) - 2*y(2) + y(3);
        if abs(den) > 0
            d = 0.5*(y(1)-y(3))/den;
            d = max(min(d,1),-1);
        end
    end
    tau(n) = (nlist(i)+d)*step;
end

if nargout > 1
    info.grid_step = step;
    info.halfrange  = Nh*step;
end
end
