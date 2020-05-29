Tini = 0;   # Tiempo inicial
Tfin = 2*pi;  # Tiempo final
n = 100;    # Instantes de tiempo donde se evaluara la solucion

function xdot = myf(x, t) 
  
  # x = (v1_x, v1_y, v2_x, v2_y, r1_x, r1_y, r2_x, r2_y)
  # Notar que aunque la variable t es  llamada no es usada
  
  # j cuenta la cantidad de veces que se llamo a la funcion
  global j; 
  j++;

  # G = 1
  m1 = 4; 
  m2 = 4;
  
  abajo = sqrt( (x(7) - x(5))^2 + (x(8) - x(6))^2 )**3;
  
  # Creo el vector para poner el return de la funcion
  xdot = zeros(8,1);
  
  # Los valores de velocidad del vector x son los nuevos valores de posicion del vector xdot
  xdot(5:8) = x(1:4);
  
  # Calculo los valores de velocidad en x y en y del vector xdot
  xdot(1) = m2*(x(7) - x(5))/abajo;
  xdot(2) = m2*(x(8) - x(6))/abajo;
  xdot(3) = m1*(x(5) - x(7))/abajo;
  xdot(4) = m1*(x(6) - x(8))/abajo;
endfunction

function u = Theta(f, theta, u0, t)
  
  # Aplica el metodo theta para resolver un sistema de ec diferenciales de primer orden.
  #
  # La idea es que si theta es cero haga el metodo EE y deje de calcular
  # y en el caso que sea distinto de cero haga el metodo CN y EI 
  # usando como semilla lo obtenido en el metodo EE
  # 
  # f     :   funcion f nombrada en la teoría
  # theta :   parametro del metodo
  # u0    :   conjunto de las condiciones iniciales (fila)
  # t     :   vector de tiempo en el que se quiere obtener la solucion (columna)
  
  n = size(t)(2);   # Cantidad de puntos en los que se requiere la solucion
  m = size(u0)(2);  # Cantidad de funciones incognitas
  
  # Matriz solucion: en cada fila tengo el valor de las incognitas,
  #                  y las columnas me marcan el paso del tiempo
  u = zeros(n, m);
  
  u(1, :) = u0;       # Primer punto : condicion inicial
  
  # Recorremos los siguientes puntos
  for i = 2:n
    #calculo paso de tiempo
    h = t(i) - t(i-1);
    
    # Euler explicito
    u(i, :) = u(i-1, :) + h*feval(f, u(i-1, :)', t(i-1))'; # Agrega una fila
    
    # Si theta es cero, entonces se termina la iteracion
    # Si theta no es cero, tengo un metodo implicito. Itero usando como semilla
    # el valor de Euler explicito.
    
    if theta == 0
      continue
    elseif theta == 0.5
      # Crank-Nicolson
      u(i, : ) = u(i-1, : ) + h*(0.5*feval(f, u(i-1, : )', t(i-1))' + 0.5*feval(f, u(i, : )', t(i))');
      continue
    endif
    
    # Euler implicito
    u(i, : ) = u(i-1, : ) + h * theta * feval(f, u(i, : )', t(i))'; 
    
  endfor
endfunction

function e = EMec(u, t)
  
  #calcula la energia mecanica del sistema de masas m1 y m2 del sistema.
  #se computa la energia por unidad de masa (dado que las masas son iguales)
  #dado el resultado de la resolucion del sistema de EDOs en la matriz u, computa
  #la energia para cada t
  
  n = size(t)(2); # Cantidad de tiempo igual a n
  
  e = zeros(n,1);
  
  for i = 1:n
    k = 2*(u(i, 1)^2 + u(i, 2)^2 + u(i, 3)^2 + u(i, 4)^2);     # Cinetica   ( 1/2 M v^2)
    p = 16 / norm(u(i, 5:6) - u(i, 7:8));                      # Potencial  ( G M1 M2 / |r1-r2| )
    e(i, 1) = k - p;               # Mecanica = Cinetica + Potencial  ( La potencial es negativa)
  endfor
  
endfunction

## feval : https://octave.sourceforge.io/octave/function/feval.html   