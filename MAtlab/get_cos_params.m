function [R, phi] =  get_cos_params(x, samples)
    N = len(samples);
    % x = np.linspace(0, 2*np.pi, N, endpoint=False)
    template = np.exp(1j * x);
    corr = 2 / N * template@samples;
    R = abs(corr);
    phi = (log(corr).imag)/(2*np.pi);
    % return R, phi