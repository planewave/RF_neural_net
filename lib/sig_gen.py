import numpy as np
# from scipy.stats import norm
import numba


def rrcosdesign(beta, span, sps):
    """
    returns the coefficients, that correspond to a square-root raised 
    cosine FIR filter with rolloff factor specified by beta. The filter 
    is truncated to span symbols, and each symbol period contains sps 
    samples. The order of the filter, sps*span, must be even. 
    modified from:
    https://github.com/veeresht/CommPy/blob/master/commpy/filters.py
    """
    n = sps*span
    rrc = np.zeros(n, dtype=float)

    for x in np.arange(n):
        t = (x-n/2)/sps
        if t == 0.0:
            rrc[x] = 1.0 - beta + (4*beta/np.pi)
        elif beta != 0 and t == 1/(4*beta):
            rrc[x] = ((beta/np.sqrt(2))*(((1+2/np.pi)*(np.sin(np.pi/
                     (4*beta)))) + ((1-2/np.pi)*(np.cos(np.pi/(4*beta))))))
        elif beta != 0 and t == -1/(4*beta):
            rrc[x] = ((beta/np.sqrt(2))*(((1+2/np.pi)*(np.sin(np.pi/
                     (4*beta)))) + ((1-2/np.pi)*(np.cos(np.pi/(4*beta))))))
        else:
            rrc[x] = ((np.sin(np.pi*t*(1-beta))+4*beta*(t)*np.cos(np.pi*
                     t*(1+beta)))/(np.pi*t*(1-(4*beta*t)*(4*beta*t))))
    return rrc


@numba.jit(nopython=True)
def upsample(x, sps, zeros=True):
    """
    increase sample rate by integer factor
    y = upsample(x,n) increases the sample rate of x by
    inserting n – 1 zeros or same sample between samples.
    input is 1D numpy array
    """
    zo = np.zeros((x.size, sps), dtype=x.dtype)
    if zeros:
        zo[:, 0] += x
    else:
        zo += np.expand_dims(x, 1)
    return zo.flatten()


@numba.jit(nopython=True)
def norm_cdf(x):
    """http://web2.uwindsor.ca/math/hlynka/zogheibhlynka.pdf"""
    return 0.5+(x+x**3/3+x**5/15+x**7/105+x**9/945)*np.exp(-x**2/2)/np.sqrt(2*np.pi)


@numba.jit(nopython=True)
def gfsk_mod(msg, sps, bt, mi):
    """
    GFSK modulator
    msg: binary message (0 or 1)
    sps: sample per symbel 
    bt: bandwidth time product
    mi: modulation index, phase change in one bit.

    the span of Gaussian filter is set to one
    """

    msg = msg*2-1.0
    freq = (np.expand_dims(msg, 1)@np.ones((1, sps))).flatten()
    t = np.arange((-0.5+1/sps/2), 0.5, 1/sps)
    shape = norm_cdf(2*np.pi*bt*(t+0.5)/np.sqrt(np.log(2)))- \
        norm_cdf(2*np.pi*bt*(t-0.5)/np.sqrt(np.log(2)))
    shape = shape/shape.sum()
    shape = (np.ones((msg.size, 1))@np.expand_dims(shape, 0)).flatten()
    freq = freq*shape
    phase = np.zeros_like(freq)+2*np.pi*np.random.random()  # random init phase
    for idx in range(freq.size-1):
        phase[idx+1] = phase[idx] + mi*np.pi*freq[idx]
    return np.exp(1j*phase)


def psk_mod(msg, M, beta, span, sps):
    """
    PSK modulator
    msg: message (0 to M-1), np.random.randint(0, M, symb)
    M: modulation order
    beta: rolloff factor of RRC filter
    span: span of RRC filter
    sps: sample per symbel
    """
    rrc = rrcosdesign(beta, span, sps)
    sig_mod = np.exp(1j*(np.pi/M+msg*(2*np.pi/M)))
    sig_up = upsample(sig_mod, sps)
    sig_pulse = np.convolve(sig_up, rrc)
    return sig_pulse[int(sps*span/2):int(1-sps*span/2)]

if __name__ == '__main__':
    pass
