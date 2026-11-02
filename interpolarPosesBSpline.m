function [Ttraj, t, info] = interpolarPosesBSpline(Tw, peso_giro, vChapa, d_acel, sigW, dt)
% Trayectoria suave por una secuencia de poses 4x4xN (waypoints).
%   Posición:     B-spline cúbica interpolante en función de s
%   Orientación:  cuaterniones con SQUAD, con el mismo s
%   Tiempo:       ley trapezoidal s(t): arranca, va a velocidad constante y frena (una sola vez en toda la pasada)
%
% Entradas
%   Tw         4x4xN   poses de paso (N >= 4). La primera y la última son la entrada y la salida, fuera de la pieza
%   peso_giro  [m/rad] cuánto pesa un giro frente a una traslación en la distancia generalizada
%   vChapa     [m/s]   velocidad que se quiere sobre la chapa (constante entre la aceleración y el frenado)
%   d_acel     [m]     distancia sobre la que se acelera (y se frena), medida sobre la chapa
%   sigW       Nx1     avance sobre la chapa en cada waypoint [m] (en la entrada y la salida, el largo recto de esos tramos)
%   dt         [s]     paso de muestreo
% Salidas
%   Ttraj      4x4xM   poses muestreadas cada dt
%   t          Mx1     tiempos
%   info       struct  sw (parámetro s de cada waypoint), s y sg (parámetro s y avance sobre la chapa de cada muestra), L, T, errPos, errRot

N = size(Tw,3);
if N < 4
    error('Hacen falta al menos 4 waypoints para la B-spline cúbica.');
end

%% Posiciones y cuaterniones de los waypoints
p = squeeze(Tw(1:3,4,:));                      % 3xN
Q = zeros(N,4);
for i = 1:N
    Q(i,:) = rot2quat(Tw(1:3,1:3,i));
end
for i = 2:N
    if dot(Q(i-1,:), Q(i,:)) < 0               % q y -q son la misma orientación: se toma el camino corto
        Q(i,:) = -Q(i,:);
    end
end

%% Parámetro s: distancia generalizada acumulada (traslación + giro)
d_pos  = vecnorm(diff(p,1,2), 2, 1)';          % (N-1)x1
d_giro = zeros(N-1,1);
for i = 1:N-1
    d_giro(i) = 2*acos(min(abs(dot(Q(i,:), Q(i+1,:))), 1));
end
d = sqrt(d_pos.^2 + (peso_giro*d_giro).^2);
if any(d < 1e-9)
    error('Hay dos waypoints repetidos.');
end
sw = [0; cumsum(d)];
L  = sw(end);

%% B-spline cúbica para la posición (pasa exacto por los waypoints)
pp = spline(sw', p);                           % spline cúbico con nudos en los s_i: es una B-spline de grado 3 interpolante

%% SQUAD: cuaterniones auxiliares
S = Q;
for i = 2:N-1
    l1 = qlog(qmul(qconj(Q(i,:)), Q(i+1,:)));
    l0 = qlog(qmul(qconj(Q(i,:)), Q(i-1,:)));
    S(i,:) = qmul(Q(i,:), qexp(-(l1 + l0)/4));
end

%% Ley de tiempo trapezoidal, sobre el avance en la chapa
% El TCP no avanza igual que la chapa (el capó es curvo y el TCP va 46 cm más arriba), así que la ley se arma con
% el avance sobre la chapa sg y después se pasa a s con una interpolación monótona. Así la pintura cae a velocidad constante.
Ls = sigW(end);
ta = 2*d_acel / vChapa;
a  = vChapa/ta;
T  = Ls/vChapa + ta;
t  = (0:dt:T)';
if T - t(end) > 1e-9
    t(end+1,1) = T;
end
sg = vChapa*(t - ta/2);                        % tramo de velocidad constante
f = t < ta;      sg(f) = 0.5*a*t(f).^2;        % arranque
f = t > T - ta;  sg(f) = Ls - 0.5*a*(T - t(f)).^2;   % frenado
sg = min(max(sg, 0), Ls);
s  = pchip(sigW, sw, sg);                      % parámetro s de la B-spline en cada instante

%% Poses muestreadas
M = numel(t);
P = ppval(pp, s')';                            % Mx3
Ttraj = repmat(eye(4), [1 1 M]);
for j = 1:M
    Ttraj(1:3,4,j)   = P(j,:)';
    Ttraj(1:3,1:3,j) = quat2rot(squadInterp(Q, S, sw, s(j)));
end

%% Controles
errPos = max(vecnorm(ppval(pp, sw') - p, 2, 1));
errRot = 0;
for i = 1:N
    Rb = quat2rot(squadInterp(Q, S, sw, sw(i)));
    errRot = max(errRot, acos(min(max((trace(Tw(1:3,1:3,i)'*Rb) - 1)/2, -1), 1)));
end
info.sw = sw;  info.s = s;  info.sg = sg;  info.L = L;  info.T = T;
info.errPos = errPos;  info.errRot = errRot;
end


%% ---- funciones auxiliares (cuaterniones [w x y z]) ----
function q = squadInterp(Q, S, sw, s)
i = find(sw <= s, 1, 'last');
i = min(max(i, 1), numel(sw)-1);
h = (s - sw(i)) / (sw(i+1) - sw(i));
a = slerp(Q(i,:), Q(i+1,:), h);
b = slerp(S(i,:), S(i+1,:), h);
q = slerp(a, b, 2*h*(1-h));
end

function q = slerp(a, b, h)
c = dot(a, b);
if c < 0
    b = -b;  c = -c;
end
if c > 0.9995
    q = a + h*(b - a);
else
    th = acos(c);
    q = (sin((1-h)*th)*a + sin(h*th)*b) / sin(th);
end
q = q / norm(q);
end

function q = rot2quat(R)
t = trace(R);
if t > 0
    s = 2*sqrt(t + 1);
    q = [s/4, (R(3,2)-R(2,3))/s, (R(1,3)-R(3,1))/s, (R(2,1)-R(1,2))/s];
else
    [~, i] = max(diag(R));
    j = mod(i,3) + 1;  k = mod(i+1,3) + 1;
    s = 2*sqrt(R(i,i) - R(j,j) - R(k,k) + 1);
    q = zeros(1,4);
    q(1)   = (R(k,j) - R(j,k)) / s;
    q(1+i) = s/4;
    q(1+j) = (R(j,i) + R(i,j)) / s;
    q(1+k) = (R(k,i) + R(i,k)) / s;
end
q = q / norm(q);
end

function R = quat2rot(q)
w = q(1);  x = q(2);  y = q(3);  z = q(4);
R = [1-2*(y^2+z^2),  2*(x*y-z*w),    2*(x*z+y*w);
     2*(x*y+z*w),    1-2*(x^2+z^2),  2*(y*z-x*w);
     2*(x*z-y*w),    2*(y*z+x*w),    1-2*(x^2+y^2)];
end

function r = qmul(a, b)
r = [a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4), ...
     a(1)*b(2) + a(2)*b(1) + a(3)*b(4) - a(4)*b(3), ...
     a(1)*b(3) - a(2)*b(4) + a(3)*b(1) + a(4)*b(2), ...
     a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + a(4)*b(1)];
end

function r = qconj(q)
r = [q(1), -q(2:4)];
end

function v = qlog(q)
u = q(2:4);  n = norm(u);
if n < 1e-12
    v = zeros(1,3);
else
    v = atan2(n, q(1)) * u / n;
end
end

function q = qexp(v)
n = norm(v);
if n < 1e-12
    q = [1 0 0 0];
else
    q = [cos(n), sin(n)*v/n];
end
end
