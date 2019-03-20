% close all
clear
clc
close all

bt = 0.3; % Bandwidth Time Product
mi = 0.5; % Modulation Index
sps = 10;

msg = [1,0,1,0,1,0,1,1, 1, 1]';
% msg = [0;1;0];
gfskMod = comm.CPMModulator('ModulationOrder', 2, ...
                            'FrequencyPulse', 'Gaussian', ... 
                            'BandwidthTimeProduct', bt, ...
                            'ModulationIndex', mi, ...
                            'BitInput', true, ...
                            'PulseLength', 1, ...
                            'SamplesPerSymbol', sps);
x = gfskMod(msg);
fq = diff(unwrap(angle(x)));
figure
plot(fq);
%%

msg_up = upsample(msg*2-1, sps);

% symetric sampling for even sps
t = (-0.5+1/sps/2):1/sps:0.5; 

% complementary cumulative distribution function
ccdf = @(x) 1-cdf(makedist('Normal'), x); 

% pulse shape of GFSK
g = ccdf(2*pi*bt*(t-1/2)/sqrt(log(2)))- ...
    ccdf(2*pi*bt*(t+1/2)/sqrt(log(2)));
g = g/sum(g);

% pulse shaping of frequency
freq = conv(msg_up, g);

% phase as integral of frequency
phase = zeros(size(freq));
for idx = 1:(size(freq)-1)
    phase(idx+1) = phase(idx) + mi*pi*freq(idx);
end
    
x = exp(1j*phase);
fq = diff(unwrap(angle(x)));
figure
plot(fq);