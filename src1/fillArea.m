function fillArea(xvalues, yvalues, color)
    x = [xvalues fliplr(xvalues)];
    y = [yvalues(:,1); flipud(yvalues(:,2))];
    fill(x, y, color, 'edgecolor', 'none', 'facealpha', 0.3,'HandleVisibility','off')
end

