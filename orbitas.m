source funciones.m
global j;

lsode_options ("integration method", "adams"); # Setear Lsode

# v1 = (0, 1), v2 =(0, -1), r1 = (1, 0), r2 = (-1, 0)
u0 = [0 1 0 -1 1 0 -1 0]; 

# Guarda el error en cada paso de tiempo para los 4 metodos
errores = zeros(n, 5);

# Guarda el tiempo computacional y el numero de evaluaciones para los 4 metodos
CostComp = zeros(2, 4);

t = linspace(Tini,Tfin,n);
errores(:,1) = t;

# j cuenta el numero de veces que se llama a la funcion f

## Euler Explicito  

j = 0; 
t0 = time();

u = Theta("myf", 0, u0, t);
CostComp(1, 1) = time()-t0;
CostComp(2, 1) = j;

for k = 2:n
  errores(k, 2) = abs(sqrt((u(k, 5))^2+u(k, 6)^2) - 1);#El valor exacto tiene una norma de 1 por ser trayectoria circular
endfor

xEE = u(:, 5);
yEE = u(:, 6);

## Crank-Nicolson

j = 0;
t0 = time();

u = Theta("myf", 0.5, u0, t);
CostComp(1, 2) = time()-t0;
CostComp(2, 2) = j;

for k = 2:n
  errores(k, 3) = abs(sqrt((u(k, 5))^2+u(k, 6)^2) - 1);#El valor exacto tiene una norma de 1 por ser trayectoria circular
endfor

xCN = u(:, 5);
yCN = u(:, 6);

## Euler Implicito

j = 0;
t0 = time();

u = Theta("myf", 1, u0, t);
CostComp(1, 3) = time()-t0;
CostComp(2, 3) = j;

for k = 2:n
  errores(k, 4) = abs(sqrt((u(k, 5))^2+u(k, 6)^2) - 1);#El valor exacto tiene una norma de 1 por ser trayectoria circular
endfor

xEI = u(:, 5);
yEI = u(:, 6);

#xEI = xEI(1:50);
#yEI = yEI(1:50);

## Lsode

j = 0;
t0 = time();

u = lsode("myf", u0, t);
CostComp(1, 4) = time()-t0;
CostComp(2, 4) = j;

for k =2:n
  errores(k, 5) = abs(sqrt((u(k, 5))^2+u(k, 6)^2) - 1);#El valor exacto tiene una norma de 1 por ser trayectoria circular
endfor

xLS = u(:, 5);
yLS = u(:, 6);

## Plot las trajectoriias

graf = plot(xEE, yEE, "-sb;EE;", xCN, yCN, "-xg;CN;", xEI, yEI, "-^r;EI;", xLS, yLS, "-.k;LS;");
xlabel("Movimiento en X");
ylabel("Movimiento en Y");
xlim([-1.7,1.7]);
ylim([-1.7,1.7]);
grid;
print "trayectorias.png";

## Guardo todo lo que hice

save "-ascii" "error.txt" errores;

costoEE = CostComp(:, 1);
costoCN = CostComp(:, 2);
costoEI = CostComp(:, 3);
costoLS = CostComp(:, 4);

save "-ascii" "CostoEE.txt" costoEE;
save "-ascii" "CostoCN.txt" costoCN;
save "-ascii" "CostoEI.txt" costoEI;
save "-ascii" "CostoLS.txt" costoLS;

## lsode    : https://octave.sourceforge.io/octave/function/lsode.html
## strjoin  : https://octave.sourceforge.io/octave/function/strjoin.html
## csvwrite : https://octave.sourceforge.io/octave/function/csvwrite.html