fig = openfig('MassFlowRate_Gulf_of_America.fig');
ax = findall(fig, 'type', 'axes');
for i = 1:length(ax)
    lines = findall(ax(i), 'type', 'line');
    for j = 1:length(lines)
        xData = get(lines(j), 'XData');
        yData = get(lines(j), 'YData');
        disp(['Line ', num2str(j), ':']);
        disp('X:'); disp(xData);
        disp('Y:'); disp(yData);
    end
end
save('GOA.mat', 'xData', 'yData');
