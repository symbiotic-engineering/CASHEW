load('JDF.mat')
figure
line(xData,yData,linewidth=2,displayname='Juan de Fuca')
load('GOA.mat')
line(xData,yData,linewidth=2,displayname='Gulf of America')
xlim([50,80])
ylim([0,40])
ylim([0,35])
legend(location='best')
xlabel('time [s]')
ylabel('Mass Flow Rate of CO2 [kg/s]')
xlabel('Time [s]')
load('GOA.mat')
y_avg = trapz(xData, yData) / (xData(end) - xData(1));
y_avg = trapz(xData, yData) / (xData(end) - xData(1));
load('JDF.mat')
y_avg = trapz(xData, yData) / (xData(end) - xData(1));
clear
figfix('Print1',12)