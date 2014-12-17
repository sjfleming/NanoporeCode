% Stephen Fleming
% July 17, 2013

% This script imports a specified data file and calls the function
% 'noiseplot.m' to make a plot of the noise.  Successive usages of this
% script will create a composite plot that can contain legend entries.

% Calls 'noiseplot.m' and 'importdata.m', which calls 'abfload.m' and/or
% 'abfloadQuiet.m'

% Uses .abf files.  I made these programs for data taken in gap-free mode.

%% Definitions

location = '/home/fleming/Documents/Physics/Golovchenko/Data/Axopatch/20130725/'; % file folder

% Create cell array of experiment file names without typing them all out
experimentsInOrder = [6 8 5 9 7]; % fill this in
experimentsInOrderS = num2str(experimentsInOrder');
names = cellstr(experimentsInOrderS);
names = strcat('2013_07_25_000',names,'.abf');

% Create cell array of legend info
volts = -200:100:200; % fill this in corresponding to order of listed experiments
%volts = [0 -200 200 -100 100];
voltsS = num2str(volts');
legends = cellstr(voltsS);
legends = strcat(legends,' mV');

figure(6)
cmap = colormap(jet(length(names)));

%% Loop through each file, load, and plot noise
for i = 1:length(names)

%% Import data
clear data
file = [location names{i}]; % putting whole file name together
startTime = 0;
endTime = 20;
[data,si] = importdata(file,startTime,endTime,1);

%% Plot the noise

[fAvg, powerAvg] = noiseplot(data,legends{i},cmap(i,:));

%% Save the power at 5kHz as a function of applied voltage

ind = find(fAvg>5000,1,'first');
powerVoltage(i) = powerAvg(ind);

end

%% Create a plot of spectral power at 5kHz as a function of applied voltage
figure(5)
semilogy(volts,powerVoltage,'ko')
title('Noise amplitude at 5kHz as a function of applied voltage')
xlabel('Applied Voltage (mV)')
ylabel('Spectral Amplitude at 5kHz')
%axis([-240 240 3*10^-5 1])
