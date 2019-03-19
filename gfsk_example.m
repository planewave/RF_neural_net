% close all
clear
clc
close all

bt = 0.8; % Bandwidth Time Product
mi = 0.5; % Modulation Index
sps_gen = 10;
% sps = 5;
% offset = 4;

n_symbo = 100;
msg = repmat([1, 0], 1, n_symbo/2)';
msg = [1,0,1,0,1,0,1,1, 1, 1]';
msg = [0;1;0];
gfskMod = comm.CPMModulator('ModulationOrder', 2, ...
                            'FrequencyPulse', 'Gaussian', ... 
                            'BandwidthTimeProduct', bt, ...
                            'ModulationIndex', mi, ...
                            'BitInput', true, ...
                            'PulseLength', 1, ...
                            'SamplesPerSymbol', sps_gen);
                        
fskMod = comm.CPMModulator('ModulationOrder', 2, ...
                            'FrequencyPulse', 'Raised Cosine', ... 
                            'BitInput', true, ...
                            'SamplesPerSymbol', sps_gen);                        
                        
x = gfskMod(msg);
% x = mskmod(msg, sps_gen);

fq = diff(unwrap(angle(x)));
figure
plot(fq);
%%

