run('robot_TF.m');

% Parámetros
pos_capo       = [2.15 0 1.2];                % [m] centro en planta y altura del borde más bajo del capot (en la terna del robot)
n_pasadas      = 5;
n_via          = 13;                          % waypoints sobre la chapa de cada pasada (la B-spline pasa exacto por ellos)
d_acel         = 0.08;                        % [m] distancia para acelerar y frenar, fuera del capot
sobre          = 0.10;                        % [m] el gatillo se abre/cierra esta distancia antes/después del borde (la pintura que se pasa del capot)
entrada        = d_acel + sobre;              % [m] largo de la entrada (y de la salida) fuera del capt
ancho_abanico  = 0.35;                        % [m] ancho del abanico sobre la chapa que contiene ±2σ (RB1000i: 250~500 mm, supuesto a confirmar)
sep_pasadas    = 0.22;                        % [m] separación entre pasadas
peso_giro      = 0.10;                        % [m/rad] cuánto pesa girar frente a trasladarse al medir el avance entre poses
dt_pas         = 0.05;                        % [s] paso de muestreo de las pasadas
x_pasadas      = pos_capo(1) + ((1:n_pasadas) - (n_pasadas+1)/2) * sep_pasadas;   % pasadas centradas en el capot: distancia de cada una al robot [m]
margen         = 0;                           % borde del capot que se deja libre [m]; 0 = la pasada va de borde a borde
d_pintura      = 0.20;                        % distancia de la cara de la campana a la chapa, sobre la normal [m]
ang_cabeza     = deg2rad(60);                 % la cabeza del RB1000i está doblada 60° respecto del eje recto (datasheet)
l_cabeza       = 0.262;                       % largo de la cabeza: del codo (TCP) al centro de la cara de la campana [m] (medido en el CAD de ABB)
w0             = [cosd(45); sind(45); 0];     % hacia dónde se inclina el eje recto (a 45° entre +x y +y del mundo); se proyecta en el plano de la chapa
v_pintura      = 0.30;                        % velocidad media sobre la chapa pintando [m/s]
suav_normal    = 1;                           % suavizado de la normal entre waypoints (sigma en waypoints; 0 = sin suavizar)
t_trans        = 3;                           % duración de cada transición [s]
n_trans        = 50;                          % puntos de cada transición
q0 = deg2rad([-53.911 -6.748 36.625 -90.473 -100.514 55.000]);   % solución de la primera pose (la entrada de la pasada 1, muñeca con q5 < 0)
q_home = deg2rad([0 0 0 0 0 0]);               % reposo: q = 0, brazo estirado; el atomizador está montado girado 180° y la cabeza queda hacia abajo

% Capó en la terna del robot
capo = stlread('capo_matlab_30k_mm.stl');
P = capo.Points / 1000;                       % mm -> m
ref = [(min(P(:,1))+max(P(:,1)))/2, (min(P(:,2))+max(P(:,2)))/2, min(P(:,3))];
% se gira 180° alrededor de z (x -> -x, y -> -y) y se traslada a pos_capo
Pw = [-(P(:,1)-ref(1)), -(P(:,2)-ref(2)), P(:,3)-ref(3)] + pos_capo;

% Superficie z = f(x,y): polinomio de grado 6 ajustado por mínimos cuadrados.
% Se usa para obtener el punto y la normal en cualquier (x,y) con derivadas suaves.
grado = 6;
[I,J] = meshgrid(0:grado);
sel = (I+J) <= grado;
I = I(sel)';  J = J(sel)';                    % exponentes (vectores fila)
xc = pos_capo(1);                             % centrado en x para que el ajuste no quede mal condicionado
A  = (Pw(:,1)-xc).^I .* Pw(:,2).^J;
c  = A \ Pw(:,3);
z_f  = @(x,y) ((x(:)-xc).^I .* y(:).^J) * c;
dzdx = @(x,y) (I .* (x(:)-xc).^max(I-1,0) .* y(:).^J) * c;
dzdy = @(x,y) (J .* (x(:)-xc).^I .* y(:).^max(J-1,0)) * c;

% 1) Pasadas: waypoints sobre la chapa 
T_pas    = cell(1,n_pasadas);                 % poses de cada pasada (4x4xM)
t_pas    = cell(1,n_pasadas);                 % tiempo de cada pose dentro de su pasada (arranca en 0)
pint_pas = cell(1,n_pasadas);                 % 1 = esa pose está sobre el capot (se pinta)

for k = 1:n_pasadas
    % largo de la pasada: el ancho del capó a esa distancia, menos el margen
    cerca = abs(Pw(:,1) - x_pasadas(k)) < 0.03;
    y_min = min(Pw(cerca,2)) + margen;
    y_max = max(Pw(cerca,2)) - margen;
    y = linspace(y_min, y_max, n_via)';
    if mod(k,2) == 0
        y = flipud(y);                        % las pasadas pares van al revés (zigzag)
    end
    x = x_pasadas(k) * ones(n_via,1);
    z = z_f(x,y);

    % normal de la superficie en cada waypoint (apunta hacia arriba)
    n = [-dzdx(x,y), -dzdy(x,y), ones(n_via,1)];
    n = n ./ vecnorm(n,2,2);

    % se suaviza la normal a lo largo de la pasada (filtro gaussiano, extremos replicados): sin esto q5 (y q4, q6)
    % hacen un cambio brusco donde la normal gira rápido, en los bordes curvos del capó
    if suav_normal > 0
        r  = ceil(3*suav_normal);
        kn = exp(-(-r:r).^2 / (2*suav_normal^2));  kn = kn / sum(kn);
        np_ = [repmat(n(1,:),r,1); n; repmat(n(end,:),r,1)];
        n  = [conv(np_(:,1),kn,'valid'), conv(np_(:,2),kn,'valid'), conv(np_(:,3),kn,'valid')];
        n  = n ./ vecnorm(n,2,2);
    end

    % pose del TCP en cada waypoint
    T_via = zeros(4,4,n_via);
    for j = 1:n_via
        nj = n(j,:)';
        w  = w0 - (w0'*nj)*nj;  w = w / norm(w);          % w0 proyectado en el plano de la chapa
        % la cabeza sale a 60° del eje recto y apunta contra la normal (el spray va hacia la chapa)
        z6 = -cos(ang_cabeza)*nj + sin(ang_cabeza)*w;     % eje recto del atomizador (a 60° del spray)
        x6 = sin(ang_cabeza)*nj + cos(ang_cabeza)*w;      % la cabeza se dobla hacia -x6 (atomizador montado girado 180° alrededor de z6)
        y6 = cross(z6, x6);
        % la campana apunta justo contra la normal (spray perpendicular a la chapa) y el TCP (codo) queda l_cabeza más arriba
        T_via(:,:,j) = [x6, y6, z6, [x(j);y(j);z(j)] + (d_pintura + l_cabeza)*nj; 0 0 0 1];
    end

    % poses de entrada y salida: fuera del capó, siguiendo la dirección de la pasada y con la misma orientación del extremo
    u_ini = T_via(1:3,4,2)   - T_via(1:3,4,1);      u_ini = u_ini / norm(u_ini);
    u_fin = T_via(1:3,4,end) - T_via(1:3,4,end-1);  u_fin = u_fin / norm(u_fin);
    T_ini = T_via(:,:,1);    T_ini(1:3,4) = T_via(1:3,4,1)   - entrada*u_ini;
    T_fin = T_via(:,:,end);  T_fin(1:3,4) = T_via(1:3,4,end) + entrada*u_fin;
    Tw = cat(3, T_ini, T_via, T_fin);

    sigW = [0; entrada + [0; cumsum(vecnorm(diff([x y z]),2,2))]; 0];   % avance sobre la chapa en cada waypoint (la entrada y la salida cuentan como rectas)
    sigW(end) = sigW(end-1) + entrada;
    [T_pas{k}, t_pas{k}, info] = interpolarPosesBSpline(Tw, peso_giro, v_pintura, d_acel, sigW, dt_pas);
    pint_pas{k} = info.sg >= sigW(2) - sobre - 1e-9 & info.sg <= sigW(end-1) + sobre + 1e-9;   % gatillo abierto sobre el capó y 'sobre' antes y después
end

% todas las pasadas juntas
filas_pas = cellfun(@numel, t_pas);           % poses de cada pasada
idx_fin = cumsum(filas_pas);
idx_ini = idx_fin - filas_pas + 1;
T_poses = cat(3, T_pas{:});
N = size(T_poses,3);
pint_poses = logical(vertcat(pint_pas{:}));

% puntos sobre la chapa (para graficar y verificar): del TCP se baja (d + l) sobre la normal, que sale de x6 y z6
x6s = squeeze(T_poses(1:3,1,:))';
z6s = squeeze(T_poses(1:3,3,:))';
nn  = -cos(ang_cabeza)*z6s + sin(ang_cabeza)*x6s;
puntos = squeeze(T_poses(1:3,4,:))' - (d_pintura + l_cabeza)*nn;

% 2) Cinemática inversa
[exito, q_tray, motivo] = ikine_propia(R, T_poses, true, q0);
if ~exito
    error('ikine_propia falló: %s', motivo);
end
% control de continuidad: ningún salto grande entre poses consecutivas de la misma pasada
dq = 0;
for k = 1:n_pasadas
    dq = max(dq, max(max(abs(diff(q_tray(idx_ini(k):idx_fin(k),:))))) * 180/pi);
end
fprintf('Salto máximo entre poses consecutivas: %.1f°\n', dq);

% 3) Trayectoria completa: pasadas (B-spline) unidas con jtraj (sin pintar)
Q = q_home;  t = 0;
tipo = 0;                                      % 0 = tramo con jtraj (espacio articular), 1 = pasada con B-spline (espacio cartesiano)
t_act = 0;

for k = 1:n_pasadas
    q_pas = q_tray(idx_ini(k):idx_fin(k), :);

    % transición hasta el primer punto de la pasada (la entrada, fuera del capó)
    tramo = jtraj(Q(end,:), q_pas(1,:), n_trans);
    tt = t_act + (1:n_trans-2)' * t_trans/(n_trans-1);
    Q = [Q; tramo(2:end-1,:)];
    t = [t; tt];
    tipo = [tipo; zeros(n_trans-2,1)];
    t_act = t_act + t_trans;

    % la pasada: pinta solo mientras el TCP está sobre el capó
    tt = t_act + t_pas{k};
    Q = [Q; q_pas];
    t = [t; tt];
    tipo = [tipo; ones(numel(tt),1)];
    t_act = tt(end);
end

% vuelta al reposo
tramo = jtraj(Q(end,:), q_home, n_trans);
tt = t_act + (1:n_trans-1)' * t_trans/(n_trans-1);
Q = [Q; tramo(2:end,:)];
t = [t; tt];
tipo = [tipo; zeros(n_trans-1,1)];

M = [t, Q];

% Verificación
err = zeros(N,1);
for i = 1:N
    Tf = R.fkine(q_tray(i,:));
    Tf = Tf.T;
    err(i) = norm(Tf(1:3,4) - T_poses(1:3,4,i));
end
fprintf('Error máximo de posición del TCP (fkine vs pose pedida): %.3f mm\n', max(err)*1000);
fprintf('Filas: %d   duración total: %.1f s\n', size(M,1), t(end));
fprintf('Velocidad articular máxima [deg/s]: %s\n', mat2str(round(rad2deg(max(abs(diff(Q)./diff(t)))),1)));
fprintf('Distancia mínima de q5 a 0 (singularidad de muñeca): %.1f°\n', rad2deg(min(abs(q_tray(:,5)))));   % solo las pasadas, sin las filas de home
v_tcp = cell(1,n_pasadas);
for k = 1:n_pasadas
    p = squeeze(T_pas{k}(1:3,4,:))';
    v_tcp{k} = vecnorm(diff(p),2,2) ./ diff(t_pas{k});
end

% Espesor de pintura sobre el capó
sigma = ancho_abanico/4;
[gx,gy] = meshgrid(min(Pw(:,1)):0.02:max(Pw(:,1)), min(Pw(:,2)):0.02:max(Pw(:,2)));
kb = boundary(Pw(:,1), Pw(:,2), 0.8);          % contorno del capó en planta
dentro = inpolygon(gx, gy, Pw(kb,1), Pw(kb,2));
gx = gx(dentro);  gy = gy(dentro);
dt_m = zeros(N,1);                             % duración de cada muestra
for k = 1:n_pasadas
    dtk = diff(t_pas{k});
    dt_m(idx_ini(k):idx_fin(k)) = [dtk; dtk(end)];
end
ip = find(pint_poses);                         % muestras con el gatillo abierto
D2 = (gx(:) - puntos(ip,1)').^2 + (gy(:) - puntos(ip,2)').^2;
esp = (exp(-D2/(2*sigma^2)) / (2*pi*sigma^2)) * (dt_m(ip) * v_pintura) * sep_pasadas;
fprintf('Espesor de pintura (abanico %.2f m): %.0f%% del capó dentro de ±15%%, %.1f%% por debajo del 50%%, %.1f%% por encima del 130%%\n', ...
    ancho_abanico, 100*mean(abs(esp-1) < 0.15), 100*mean(esp < 0.5), 100*mean(esp > 1.3));

% Gráficos
figure('Name','Articulaciones','NumberTitle','off');
for j = 1:6
    subplot(3,2,j);
    plot(t, rad2deg(Q(:,j)), 'LineWidth', 1.5); hold on;
    yline(rad2deg(R.links(j).qlim(1)), 'r--');
    yline(rad2deg(R.links(j).qlim(2)), 'r--');
    title(sprintf('q%d', j)); xlabel('t [s]'); ylabel('[°]'); grid on;
end
sgtitle('Trayectoria articular completa');

figure('Name','Velocidad del TCP','NumberTitle','off'); hold on;
for k = 1:n_pasadas
    plot(t_pas{k}(2:end), v_tcp{k}, 'LineWidth', 1.5);
end
grid on; xlabel('t [s] (desde el inicio de cada pasada)'); ylabel('[m/s]');
legend(compose('pasada %d', 1:n_pasadas)); title('Velocidad del TCP en cada pasada');

% Robot con trail, tramo por tramo: jtraj (rojo) y pasadas con B-spline (azul)
corte = [1; find(diff(tipo) ~= 0) + 1; numel(tipo) + 1];
n_tramos = numel(corte) - 1;

% recorrido del TCP de cada tramo (estático, para ver de un vistazo qué hace cada función)
tcp_Q = zeros(size(Q,1), 3);
for i = 1:size(Q,1)
    Tf = R.fkine(Q(i,:)).T;
    tcp_Q(i,:) = Tf(1:3,4)';
end
figure('Name','TCP por tramos','NumberTitle','off');
trisurf(capo.ConnectivityList, Pw(:,1), Pw(:,2), Pw(:,3), 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha',0.5);
hold on; axis equal; grid on; view(-35,30);
for s = 1:n_tramos
    ii = max(corte(s)-1,1):corte(s+1)-1;
    if tipo(corte(s)) == 1, c = 'b'; else, c = 'r'; end
    plot3(tcp_Q(ii,1), tcp_Q(ii,2), tcp_Q(ii,3), c, 'LineWidth', 1.5);
end
hr = plot3(NaN,NaN,NaN,'r','LineWidth',1.5);  hb = plot3(NaN,NaN,NaN,'b','LineWidth',1.5);
legend([hr hb], {'jtraj (articular)','B-spline (cartesiano)'});
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]'); title('Recorrido del TCP por tramos');

% animación del robot con 'trail' en cada tramo, sobre el capot y las pasadas
fig_anim = figure('Name','Animación por tramos','NumberTitle','off');
trisurf(capo.ConnectivityList, Pw(:,1), Pw(:,2), Pw(:,3), 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha',0.5);
hold on; grid on;
for s = 1:n_tramos
    ii = max(corte(s)-1,1):corte(s+1)-1;
    if tipo(corte(s)) == 1, c = 'b'; else, c = 'r'; end
    plot3(tcp_Q(ii,1), tcp_Q(ii,2), tcp_Q(ii,3), c, 'LineWidth', 1.5);
end
try
    for s = 1:n_tramos
        figure(fig_anim);                             % fuerza que el robot se dibuje en esta figura
        ii = max(corte(s)-1,1):corte(s+1)-1;
        ii = unique([ii(1:2:end), ii(end)]);          % se saltea un cuadro de cada dos
        if tipo(corte(s)) == 1, estilo = 'b-'; else, estilo = 'r-'; end
        R.plot(Q(ii,:), 'workspace', [-1 3.5 -2 2 -0.5 3], 'fps', 15, 'trail', estilo);
    end
catch
    disp('Animación interrumpida (se cerró la figura); se sigue con el resto del script.');
end

% Diagnóstico: el spray sale perpendicular a la chapa? 
TR = triangulation(capo.ConnectivityList, Pw);
nF = faceNormal(TR);
nF(nF(:,3)<0,:) = -nF(nF(:,3)<0,:);           % normales hacia arriba
C  = incenter(TR);
spray_dir = cos(ang_cabeza)*z6s - sin(ang_cabeza)*x6s;    % dirección de la cabeza (por donde sale la pintura)
radio = 0.03;
ang = nan(N,1);
for i = find(pint_poses)'
    d = vecnorm(C - puntos(i,:), 2, 2);
    cerca = d < radio;
    nm = mean(nF(cerca,:), 1);  nm = nm / norm(nm);
    ang(i) = acosd(-spray_dir(i,:) * nm');
end
fprintf('Desvío del spray contra la normal de la malla: medio %.2f°, máximo %.2f°\n', mean(ang,'omitnan'), max(ang));
