clc
clear
global a
a = 6378137; 
global e2
e2 = 0.00669437999013; 
%wsp samolotu 
macierzDane = load('flight.txt'); 
phi = macierzDane(:,1)*pi()/180;  %zamiana na rad
lambda = macierzDane(:,2)*pi()/180;  %zamiana na rad
h = macierzDane(:,3); 
%wsp lotniska
phiB = phi(1);
lambdaB = lambda(1);
hB = h(1);

neu = [0, 0, 0];

for i=2:1909
    Xij = geo_neu(phi(i), lambda(i), h(i), phiB, lambdaB, hB);
    neu = [neu; Xij];
    i = i+1;
end

n = neu(:,1);
e = neu(:,2);
u = neu(:,3);

[A1, s1, z1] = neu_asz(n(1), e(1), u(1));
Asz = [A1, s1, z1];

for i=2:1909
    [Ai, si, zi] = neu_asz(n(i), e(i), u(i));
    Asz = [Asz; Ai, si, zi];
end

i = 2;
while i <= 1909
    if u(i) <= 0, break;end
    i=i+1;
end
wspolrzedne_neu_gdy_samolot_znika_za_horyzontem = neu(3,:)

% rysowanie
geoscatter(phi*180/pi,lambda*180/pi,5, 'ro')
plot3(n, e, u)



function N = liczN(fi)
    global a
    global e2
    N = a / sqrt(1 - e2*((sin(fi))^2));
end

function [x, y, z] = geo_xyz(f, l, h)
    global e2
    N = liczN(f);
    x = (N + h) * cos(f) * cos(l);
    y = (N + h) * cos(f) * sin(l);
    z = ((N * (1 - e2)) + h) * sin(f);
end

function Xij = geo_neu(f1, l1, h1, f2, l2, h2)
    [a1, a2, a3] = geo_xyz(f1, l1, h1);
    [b1, b2, b3] = geo_xyz(f2, l2, h2);
    vec_delt = [b1-a1; b2-a2; b3-a3];

    Xij = [-sin(f1) * cos(l1), -sin(l1), cos(f1) * cos(l1);
          -sin(f1) * sin(l1), cos(l1), cos(f1) * sin(l1);
          cos(f1), 0, sin(f1)];
    Xij = Xij * vec_delt;
    Xij = Xij';
end

function [A, s, z] = neu_asz(n, e, u)
    A = atan(e/n);
    s = sqrt(n*n + e*e + u*u);
    z = acos(u/s);
end
