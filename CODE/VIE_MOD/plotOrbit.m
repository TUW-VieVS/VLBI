function [ ] = plotOrbit(sourcesO, sourcesN)

rO = [sourcesO.s.x_crf, sourcesO.s.y_crf, sourcesO.s.z_crf];
rN = [sourcesN.s.x_crf, sourcesN.s.y_crf, sourcesN.s.z_crf];

rO = rO(1:168,:);
rN = rN(1:168,:);
figure()
scatter3(rO(:,1), rO(:,2), rO(:,3), '.')
hold on
scatter3(rO(1,1), rO(1,2), rO(1,3), 100, '*', 'blue')
%quiver3(0,0,0, n(1)/norm(n) *4e7, n(2)/norm(n) *4e7, n(3)/norm(n) *4e7, 'green', 'LineWidth',2); %RAAN
%plot3([0, h(1)/10000],  [0, h(2)/10000],  [0, h(3)/10000],'blue')
%quiver3(0,0,0, x(1)/norm(x) *3e7, x(2)/norm(x) *3e7, x(3)/norm(x) *3e7,'magenta', 'LineWidth',2) 
%quiver3(0,0,0, y(1)/norm(y) *2.5e7, y(2)/norm(y) *2.5e7, y(3)/norm(y) *2.5e7,'magenta', 'LineWidth',2) %y-achse
%quiver3(0,0,0, z(1)/norm(z) *4e7, z(2)/norm(z) *4e7, z(3)/norm(z) *4e7,'magenta', 'LineWidth',2) %z-Achse

scatter3(rN(:,1), rN(:,2), rN(:,3), '.', 'red')
scatter3(rN(1,1), rN(1,2), rN(1,3), 100, '*', 'red')
end