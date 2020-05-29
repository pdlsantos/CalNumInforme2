source funciones.m
global j;

lsode_options ("integration method","adams"); # Setear Lsode

# v1 = (0, 0.5), v2 =(0, -0.5), r1 = (1, 0), r2 = (-1, 0)
u0 = [0 0.5 0 -0.5 1 0 -1 0];

# Tiene un período distinto
Tfin = 2.68;

t = linspace(Tini, Tfin, n);

e = zeros(n,  5);
e(:,1) = t;

## Euler Explicito  

u = Theta("myf", 0, u0, t);
e(:,2) = EMec(u,  t);

## Crank-Nicolson

u = Theta("myf", 0.5, u0,  t);
e(:,3) = EMec(u, t);

## Euler Implicito
u = Theta("myf", 1, u0, t);
e(:,4) = EMec(u, t);

#lsode
u = lsode("myf", u0, t);
e(:,5) = EMec(u, t);

# Comprueba que el periodo esta aproximadamente bien
graf = plot(u(:,5), u(:,6), ".r");
xlabel("X");
ylabel("Y");
xlim([-0.2,1.2]);
ylim([-0.5,0.5]);
grid;

## guardar los valores de energia
save "-ascii" "EnergiaMec.txt" e;

## lsode : https://octave.sourceforge.io/octave/function/lsode.html