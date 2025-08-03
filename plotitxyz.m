function plotitxyz(x, y, z, xm, ym, zm, imin,Re,Rm)

figure('Name', 'Spacecraft trajectory in Moon-fixed rotating frame', ...
       'Color', [1 1 1]);

[xx, yy, zz] = sphere(128);
hold on

%... Spacecraft trajectory:
plot3(x, y, z, 'r', 'LineWidth', 2.0)

%... Moon trajectory:
plot3(xm, ym, zm, 'g', 'LineWidth', 0.5)

%... Earth:
Earth = surfl(Re * xx, Re * yy, Re * zz);
set(Earth, 'FaceAlpha', 0.5);
shading interp

%... Geocentric moon-fixed coordinate axes:
L1 = 63 * Re; 
L2 = 20 * Re; 
L3 = 29 * Re;

line([0 L1], [0 0], [0 0], 'color', 'k')
text(L1, 0, 0, 'x', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino')

line([0 0], [0 L2], [0 0], 'color', 'k')
text(0, L2, 0, 'y', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino')

line([0 0], [0 0], [0 L3], 'color', 'k')
text(0, 0, L3, 'z', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino')

%... Spacecraft at TLI
plot3(x(1), y(1), z(1), 'o', ...
      'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 3)

%... Spacecraft at closest approach
plot3(x(imin), y(imin), z(imin), 'o', ...
      'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 2)

%... Spacecraft at tf (end)
plot3(x(end), y(end), z(end), 'o', ...
      'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 3)

%... Moon at TLI
text(xm(1), ym(1), zm(1), 'Moon at TLI')
Moon1 = surfl(Rm * xx + xm(1), Rm * yy + ym(1), Rm * zz + zm(1));
set(Moon1, 'FaceAlpha', 0.99)
shading interp

%... Moon at closest approach
Moon2 = surfl(Rm * xx + xm(imin), Rm * yy + ym(imin), Rm * zz + zm(imin));
set(Moon2, 'FaceAlpha', 0.99)
shading interp

%... Moon at end of simulation
Moon3 = surfl(Rm * xx + xm(end), Rm * yy + ym(end), Rm * zz + zm(end));
set(Moon3, 'FaceAlpha', 0.99)
shading interp

axis image
axis vis3d
axis off
view([1, 1, 1])
end
