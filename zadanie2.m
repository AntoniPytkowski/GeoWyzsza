% https://stellarium-web.org/skysource/GaiaDR21114998150271095040?fov=5.7818&date=2022-01-18T16:34:21Z&lat=52.25&lng=21.00&elev=0
clear
clc

% dane gwiazdy Gaia DR2
% 06h 40m44.0s   +73°40'37.7"
re = 6 + 40/60 + 44/3600;
de = 73 + 40/60 + 37.7/3600;

% lokalizacja obserwatora (las bielanski)
fi = 52.293455;
la = 20.972900;

tab = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24];
% tab to kolejne godziny w dobie

for i = 1:24
    tab(i) = mod(katgodz(2022, 1, 18, tab(i), la, re), 360);
end
% tab to kat godzinny (< 360) dla kolejnych godzin

cosZ = (sind(fi)*sind(de) + cosd(fi)*cosd(de).*cosd(tab));
Z = acosd(cosZ);
h = 90 - Z;

dy = -cosd(de).*sind(tab);
dx = cosd(fi).*sind(de) - sind(fi).*cosd(de).*cosd(tab);
az = mod(atan2d(dy, dx), 360);

x = sind(Z).*cosd(az);
y = sind(Z).*sind(az);
z = cosd(Z);

% wizualizacja 
[ax, ay, az] = sphere(24);
ax = ax(13:end,:);
ay = ay(13:end,:);
az = az(13:end,:);
surf(ax,ay,az,'FaceColor','blue','FaceAlpha',0.5)
axis equal
hold on

scatter3(x,y,z, 15, 'yellow', '*')


function [t] = katgodz(y,d,m,h,lambda,alfa) 
jd = juliandate(datetime(y,m,d)); %dni 
g = GMST(jd); %stopnie 
UT1 = h*1.002737909350795; %godziny 
%obliczenie czasu gwiazdowego(w stopniach) 
S = UT1*15 + lambda + g; 
%obliczenie kąta godzinnego(w stopniach) 
t = S - alfa*15; 
end

function g = GMST(jd)
    T = (jd - 2451545) / 36525;
    g = 280.46061837 + 360.98564736629  * (jd - 2451545.0) + 0.000387933*T.^2 - T.^3/38710000;
    g = mod(g, 360);
end