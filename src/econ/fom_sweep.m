clear;clc;close all

eta = linspace(0.4,1); % efficiency
lcoe = linspace(0.01,.5)/1e3; % levelized (per year) dollar per Wh

C_ann_other = [25 50 100]*1e5; % levelized (per year) dollars
Wdot = 1e6; % 1 MW

figure
for i=1:length(C_ann_other)
    [ETA,LCOE] = meshgrid(eta,lcoe);
    FOM = 1 ./ ETA .* (LCOE + C_ann_other(i)/(Wdot*8766));
    
    sp = subplot(1,length(C_ann_other),i);
    contourf(ETA,1e3*LCOE,FOM,'DisplayName','Figure of Merit')
    
    caxis([0.3 4]*1e-3)
    xlabel('Pump Efficiency')
    ylabel('Levelized cost of energy ($/kWh)')
    title(['C_{ann,other}=' sprintf('%.1e', C_ann_other(i))])

    triangle = native2unicode([0xF0 0x9F 0x9F 0x80],'utf8');
    text(.5, .25,triangle,'Color',[0 0 0],'FontSize', 25)
    hold on
    plot(1,.5,'ro','MarkerFaceColor','r','DisplayName','Wave')
    plot(.5,.03,'o','Color',[1 1 0],'MarkerFaceColor',[1 1 0],'DisplayName','Onshore Wind or Solar')
    plot(NaN,NaN,'ko','DisplayName','Floating Offshore Wind')
    if i==length(C_ann_other)
        legend
        colorbar
    end

    position_plot(sp)

end
improvePlot

% from https://www.mathworks.com/matlabcentral/answers/2080611-how-to-make-widths-of-all-subplot-and-colorbars-same
function position_plot(ax)
    set(ax,'Units','normalized');
    set(ax,'Position',[ax.Position(1), ax.Position(2), 0.2, ax.Position(4)]); % You can change 0.4 (the width) to meet your needs
end