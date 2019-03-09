import numpy as np

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


def gaussdesign(bt, span, sps):
    """
    Gaussian FIR Pulse-Shaping Filter Design
    converted from Matlab
    
    pulse-shaping filter. BT is the 3-dB bandwidth-symbol time product,
    where B is the one-sided bandwidth in hertz and T is the symbol time in
    seconds. The filter is truncated to SPAN symbols and each symbol is
    represented by SPS samples.
    the filter order which is SPS*SPAN, must be even
    """

    filtLen = sps*span+1
    t = np.linspace(-span/2,span/2,filtLen)
    alpha = np.sqrt(np.log(2)/2)/(bt)
    h = (np.sqrt(np.pi)/alpha)*np.exp(-(t*np.pi/alpha)**2)
    h = h/np.sum(h)
    return h


def upsample(x, sps, zeros=True):
    """
    increase sample rate by integer factor
    y = upsample(x,n) increases the sample rate of x by 
    inserting n – 1 zeros or same sample between samples.
    input is 1D numpy array
    
    """
    zo = np.zeros((x.size, sps), dtype=x.dtype)
    if zeros:
        zo[:,0] += x
    else:
        zo += np.expand_dims(x, 1)
            
    return zo.flatten()

def psk_gen(symb=64, M=4, beta=0.4, span=4, sps=16):
    """
    generate baseband PSK signal
    """
    rrc = rrcosdesign(beta, span, sps)
    msg = np.random.randint(0, M, symb) 
    sig_mod = np.exp(1j*(np.pi/M+msg*(2*np.pi/M)))
    sig_up = upsample(sig_mod, sps)
    sig_pulse = np.convolve(sig_up, rrc)
    return sig_pulse[int(sps*span/2):int(1-sps*span/2)]


def gfsk_gen(symb=64, bt=0.6, mi=0.5, span=4, sps=16):
    """
    generate baseband GFSK signal
    bt: bandwidth time product
    mi: modulation index, phase change in one bit.
        when mi = 0.5 is MSK. 0.5 pi change per symbol
    """
    gaus = gaussdesign(bt, span, sps)
    msg = np.random.randint(0, 2, symb)
    freq = msg*2-1.0
    freq = upsample(freq, sps, zeros=False)
    freq_gaus = np.convolve(freq, gaus)
    phase = np.zeros_like(freq_gaus)
    for idx in range(freq_gaus.size-1):
        phase[idx+1] = phase[idx] + np.pi/sps*mi*freq_gaus[idx]
    sig = np.exp(1j*phase)
    return sig[int(sps*span/2):int(-sps*span/2)]

if __name__ == '__main__':
    pass
