 T1 = readtable('acceleration_250Hz.csv'); %PPG data
 T2 = readtable('ECG_120Hz.csv'); %Combined ACC and ECG (was pt/mt
% cd(thisDir)

T1Arr = table2array(T1); %Break into arrays
T2Arr = table2array(T2);

plot(T1Arr(:,1),T1Arr(:,2))
figure

plot(T2Arr(:,1),T2Arr(:,2))