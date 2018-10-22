function [phi,lam,H] = kart2ell(P,a,b)

ee = (a^2-b^2)/a^2;

'Entfernung des Punktes von der z-Achse';
r1 = sqrt(P(1)^2+P(2)^2);

'Näherungswert für beta';
beta = atan(P(3)/r1);

'Berechnung der Hilfsgrößen';
m = b*P(3)/(a*r1);
n = (a^2-b^2)/(a*r1);

'Bestimmung der reduzierten Breite beta';
for i=1:5
    beta = atan(m+n*sin(beta));
end

'Bestimmung der geodätischen Breite phi';
phi = atan(a/b*tan(beta));

'Bestimmung der Länge lam';
lam = atan2(P(2),P(1));

'Querkrümmungsradius N';
N = a/sqrt(1-ee*sin(phi)^2);

'ellipsoidische Höhe H';
H = P(3)/sin(phi)-(1-ee)*N;
