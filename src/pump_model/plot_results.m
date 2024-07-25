addpath ../sea-lab-utils/plotutilities/
get_colors();       % get color blind friendly colors

% Extract Data
F =  out.simout.Force_on_WEC.Data;
mdot = out.simout.Sequestered_CO2.Data;
time = out.tout;

% only plot once system has settled
ii = 1;
while time(ii)<5
    ii = ii+1;
end
time = time(ii:end) - time(ii);
F = F(ii:end)/1e6;      % [MN]
mdot = mdot(ii:end);

figure('Name','simscape results')
subplot(2,1,1)
line(time,mdot,'color',blue,'linewidth',2)
xticklabels([])
ylabel({'Sequestration [kg/s]'})
ylim([0,300])
figfix('Print1',7);
subplot(2,1,2)
line(time,F,'color',blue,'linewidth',2)
xlabel('time [s]')
ylabel('Force on WEC [MN]')
figfix('Print1',7);
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) pos(2)+0.02 pos(3) pos(4)]);