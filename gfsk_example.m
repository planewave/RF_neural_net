% close all
clear
clc


bt = 0.3; % Bandwidth Time Product
mi = 0.5; % Modulation Index
sps_gen = 10;
% sps = 5;
% offset = 4;

n_symbo = 100;
msg = repmat([1, 0], 1, n_symbo/2)';
msg = [1,0,1,0,1,0,1,1, 1, 1]';
gfskMod = comm.CPMModulator('ModulationOrder', 2, ...
                            'FrequencyPulse', 'Gaussian', ... 
                            'BandwidthTimeProduct', bt, ...
                            'ModulationIndex', mi, ...
                            'BitInput', true, ...
                            'SamplesPerSymbol', sps_gen);
                        
fskMod = comm.CPMModulator('ModulationOrder', 2, ...
                            'FrequencyPulse', 'Raised Cosine', ... 
                            'BitInput', true, ...
                            'SamplesPerSymbol', sps_gen);                        
                        
x = fskMod(msg);
% x = mskmod(msg, sps_gen);

fq = diff(unwrap(angle(x)));
figure
plot(fq);

% plot(unwrap(angle(x))/pi)
grid
x = awgn(x, 23);