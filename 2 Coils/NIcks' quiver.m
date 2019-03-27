    %Plot the three solenoids
    plot3(loop(1, :), loop(2, :), loop(3, :), 'LineWidth', 1.5, 'Color', 'r');
    hold on;
    plot3(loop2(1, :), loop2(2, :), loop2(3, :), 'LineWidth', 1.5,'Color', 'b');
    plot3(loop3(1, :), loop3(2, :), loop3(3, :), 'LineWidth', 1.5,'Color', 'g');
    axis equal; %ensures axis is square to enhance clarity
    
    %Add title and labels to axis
    title({'Quiver Plot of Magnetic','Field of Helmholtz Coils'})
    xlabel('X Position (m) =>', 'FontSize', 10)
    ylabel('<= Y Position (m)', 'FontSize', 10)
    zlabel('Z Position (m) =>', 'FontSize', 10)
    
    %Create scaled unitary vectors to quiver plot
    unit_B = ones(size(B));
    for i = 1:length(B(1, :))
        unit_B(:, i) = B(:, i)/norm(B(:, i))*scale;
    end
    
    %Limits the number of vectors plotted
    coords = coordinates(:, 1:quiver_step:end);
    unit_B = unit_B(:, 1:quiver_step:end);
    
    %Plots quiver plot
    q = quiver3(coords(1, :), coords(2, :), coords(3, :), unit_B(1, :), unit_B(2, :), unit_B(3, :), 'AutoScale', 'off', 'Color', [0.75 0.75 0.75]);
    q.LineWidth = 1.5;
    grid on;
    hold off;