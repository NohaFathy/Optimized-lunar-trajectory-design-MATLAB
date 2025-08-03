function plotit_XYZ(X, Y, Z, Xm, Ym, Zm, imin, Re, Rm)
%–––––––––––––––––––––––––––––––––––––

figure('Name','Trajectories of Spacecraft (red) and Moon (green)', 'Color', 'w');


[xx, yy, zz] = sphere(128);

hold on

%... Geocentric inertial coordinate axes
L = 20 * Re;
line([0 L], [0 0], [0 0], 'color','b')
text(L, 0, 0, 'X', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino')
line([0 0], [0 L], [0 0], 'color','b')
text(0, L, 0, 'Y', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino')
line([0 0], [0 0], [0 L], 'color','b')
text(0, 0, L, 'Z', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino')

%... Earth
Earth = surfl(Re*xx, Re*yy, Re*zz);
set(Earth, 'FaceAlpha', 0.5);
shading interp

%... Spacecraft at TLI
plot3(X(1), Y(1), Z(1), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 5)

%... Spacecraft at closest approach (perilune)
plot3(X(imin), Y(imin), Z(imin), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4)

%... Spacecraft at final time
plot3(X(end), Y(end), Z(end), 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5)

%... Moon at TLI
text(Xm(1), Ym(1), Zm(1), 'Moon at TLI')
Moon1 = surfl(Rm*xx + Xm(1), Rm*yy + Ym(1), Rm*zz + Zm(1));
set(Moon1, 'FaceAlpha', 0.99)
shading interp

%... Moon at closest approach
Moon2 = surfl(Rm*xx + Xm(imin), Rm*yy + Ym(imin), Rm*zz + Zm(imin));
set(Moon2, 'FaceAlpha', 0.99)
shading interp

%... Moon at end of simulation
Moon3 = surfl(Rm*xx + Xm(end), Rm*yy + Ym(end), Rm*zz + Zm(end));
set(Moon3, 'FaceAlpha', 0.99)
shading interp

%... Spacecraft trajectory
%plot3(X, Y, Z, 'r', 'LineWidth', 1.5)

%... Moon trajectory
%plot3(Xm, Ym, Zm, 'g', 'LineWidth', 1)
% Create animated lines
spacecraftLine = animatedline('Color','r','LineWidth',1.5);
moonLine       = animatedline('Color','g','LineWidth',1.2);

% Animation loop
for k = 1:length(X)
     addpoints(spacecraftLine, X(k), Y(k), Z(k));
    addpoints(moonLine, Xm(k), Ym(k), Zm(k));
    

    h1 = plot3(X(k), Y(k), Z(k), '^', 'MarkerFaceColor', 'r', ...
               'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    
  
    t1 = text(X(k), Y(k), Z(k) + 0.1*Re, 'Spacecraft', ...
              'FontSize', 10, 'FontWeight','bold', ...
              'Color', 'r', 'HorizontalAlignment','center');

    
    h2 = plot3(Xm(k), Ym(k), Zm(k), 'o', 'MarkerFaceColor', 'g', ...
               'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    
    
    t2 = text(Xm(k), Ym(k), Zm(k) + 0.1*Re, 'Moon', ...
              'FontSize', 10, 'FontWeight','bold', ...
              'Color', 'g', 'HorizontalAlignment','center');
    
    drawnow limitrate

    % Remove current markers and labels
    if k < length(X)
        delete(h1); delete(h2);
        delete(t1); delete(t2);
    end
end
 









