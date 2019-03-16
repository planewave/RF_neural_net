close all
clear
clc


bt = 0.5; % Bandwidth Time Product
mi = 0.5; % Modulation Index
sps_gen = 50;
% sps = 5;
% offset = 4;

n_symbo = 100;
msg = repmat([1, 0], 1, n_symbo/2)';
% msg = randi([0 1],n_symbo,1);
gfskMod = comm.CPMModulator('ModulationOrder', 2, ...
                            'FrequencyPulse', 'Gaussian', ... 
                            'BandwidthTimeProduct', bt, ...
                            'ModulationIndex', mi, ...
                            'BitInput', true, ...
                            'SamplesPerSymbol', sps_gen);
x = gfskMod(msg);
% x = mskmod(msg, sps_gen);

phase = diff(unwrap(angle(x)));
plot(phase/pi);
grid
x = awgn(x, 23);