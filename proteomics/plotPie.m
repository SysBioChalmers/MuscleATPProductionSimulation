function  plotPie(results, labels)
    exMap = [67 116 160
             80 137 188
             91 155 213
             151 185 224
             190 209 234]/255;
    colormap(exMap)
    p = pie(results./sum(results), labels);
    for i = 1:length(p)
        if not(strcmp(class(p(i)), 'matlab.graphics.primitive.Text'))
            p(i).EdgeColor = 'w';
        else
            p(i).FontSize = 15;
        end
    end

end

