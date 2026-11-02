%Resolvemos la CI numerica para una sucesion de poses
%Entradas: robot, poses, q_mejor (0 o 1), q0, mascara (todas activas)
%Salidas: exito, resultado, motivo
function [exito, resultado, motivo] = ikine_propia(robot, poses, q_mejor, q0, mascara)

motivo    = '';
resultado = [];
n_joints  = robot.n;

%Si no se pasa mascara, se usan las 6 componentes de la pose
if nargin < 5 || isempty(mascara)
    mascara = [1 1 1 1 1 1];
end

%Validamos el formato correcto de la mascara
if numel(mascara) ~= 6 || ~all(mascara == 0 | mascara == 1)
    exito  = false;
    motivo = 'La máscara debe ser un vector de 6 elementos binarios (0 o 1).';
    return
end
mascara = mascara(:)';

T_poses = convertir_poses(poses);
N = size(T_poses, 3);

if N == 0
    exito  = false;
    motivo = 'La sucesión de poses está vacía.';
    return
end

%Verificamos que q0 tenga la cantidad correspondiente de elementos
if numel(q0) ~= n_joints
    exito  = false;
    motivo = sprintf('q0 debe tener %d elementos (uno por articulación).', n_joints);
    return
end
q0 = q0(:)';

if q_mejor
%Se resuelve pose por pose usando la solucion anterior como semilla inicial
%Esto garantiza continuidad en el espacio articular
    trayectoria = zeros(N, n_joints);
    q_semilla   = q0;

    for i = 1:N
        T_i = squeeze(T_poses(:,:,i));
        q_i = robot.ikine(T_i, 'q0', q_semilla, 'mask', mascara);
        %Verificamos convergencia
        if isempty(q_i) || any(isnan(q_i))
            exito     = false;
            motivo    = sprintf('ikine no convergió en la pose %d de %d.', i, N);
            resultado = trayectoria(1:i-1, :);
            return
        end
        %ikine puede caer en la otra rama de la muñeca (q4+pi, -q5, q6+pi) que da la
        %misma pose; nos quedamos con la que esté más cerca de la solución anterior
        q_i = elegir_rama_continua(robot, q_i, q_semilla);
        %Verificamos que la solucion este dentro de los limites articulares
        if ~dentro_de_limites(robot, q_i)
            exito     = false;
            motivo    = sprintf('Solución fuera del espacio articular en pose %d.', i);
            resultado = trayectoria(1:i-1, :);
            return
        end

        trayectoria(i, :) = q_i;
        q_semilla = q_i; %La solucion actual es la semilla de la siguiente pose
    end

    exito     = true;
    resultado = trayectoria;

else
    %Se explora el espacio articular con semillas aleatorias para encontrar
    %todas las soluciones distintas de la primera pose. Cada una inicia una
    %rama que se propaga por el resto de la trayectoria
    N_intentos   = 300; %Cantidad de semillas aleatorias de la primera pose
    tol_repetida = 0.05;

    T_primera    = squeeze(T_poses(:,:,1));
    sols_primera = buscar_todas_soluciones(robot, T_primera, mascara, N_intentos, tol_repetida);

    if isempty(sols_primera)
        exito  = false;
        motivo = 'No se encontró ninguna solución válida para la primera pose.';
        return
    end
    %Inicializamos una rama por cada solucion encontrada
    n_ramas       = size(sols_primera, 1);
    ramas         = cell(n_ramas, 1);
    ramas_validas = true(n_ramas, 1);

    for k = 1:n_ramas
        ramas{k}      = zeros(N, n_joints);
        ramas{k}(1,:) = sols_primera(k,:);
    end
    %Propagamos cada rama usando continuidad articular
    for i = 2:N
        T_i = squeeze(T_poses(:,:,i));

        for k = 1:n_ramas
            if ~ramas_validas(k)
                continue
            end

            q_semilla_k = ramas{k}(i-1, :);
            q_i = robot.ikine(T_i, 'q0', q_semilla_k, 'mask', mascara);
            %Si la rama falla en alguna pose se descarta
            if isempty(q_i) || any(isnan(q_i)) || ~dentro_de_limites(robot, q_i)
                ramas_validas(k) = false;
                continue
            end

            ramas{k}(i,:) = q_i;
        end
    end
    %Conservamos solo las ramas que completaron toda la trayectoria
    resultado = {};
    for k = 1:n_ramas
        if ramas_validas(k)
            resultado{end+1} = ramas{k}; %#ok<AGROW>
        end
    end

    if isempty(resultado)
        exito  = false;
        motivo = 'Ninguna rama pudo propagarse por toda la trayectoria.';
        return
    end
   %Elimiamos trayectorias que sean casi identicas
    resultado = eliminar_repetidas(resultado, tol_repetida);
    exito = true;

end
end

function q_sel = elegir_rama_continua(robot, q, q_ref)
%De las dos posturas de muñeca que dan la misma pose (la dada y la invertida)
%elige la más cercana a q_ref, trayendo cada ángulo a la vuelta que mejor entra en qlim
cand = [q; q(1:3), q(4)+pi, -q(5), q(6)+pi];
q_sel = q;  mejor = inf;
for c = 1:size(cand,1)
    qc = cand(c,:);  valido = true;
    for j = 1:robot.n
        lim = robot.links(j).qlim;
        opciones = qc(j) + 2*pi*(-2:2);
        opciones = opciones(opciones >= lim(1) & opciones <= lim(2));
        if isempty(opciones), valido = false; break, end
        [~, k] = min(abs(opciones - q_ref(j)));
        qc(j) = opciones(k);
    end
    if valido
        d = norm(qc - q_ref);
        if d < mejor, mejor = d; q_sel = qc; end
    end
end
end

function T_array = convertir_poses(poses)
%Convierte el argumento poses a matrices homogeneas
if ~isnumeric(poses)
    error('poses debe ser una matriz numérica 4x4 o 4x4xN.');
end
if ndims(poses) == 3 && size(poses,1)==4 && size(poses,2)==4
    T_array = poses;
elseif ndims(poses) == 2 && size(poses,1)==4 && size(poses,2)==4
    T_array = reshape(poses, 4, 4, 1);
else
    error('poses debe ser una matriz 4x4 (una pose) o 4x4xN (varias poses).');
end
end


function ok = dentro_de_limites(robot, q)
%Devuelve true si todos los angulos articulares estan dentro de qlim
ok = true;
for j = 1:robot.n
    if q(j) < robot.links(j).qlim(1) || q(j) > robot.links(j).qlim(2)
        ok = false;
        return
    end
end
end


function soluciones = buscar_todas_soluciones(robot, T, mascara, N_intentos, tol)
%Se usa ikine con semillas aleatorias dentro del espacio articular
n          = robot.n;
soluciones = zeros(0, n);

qlim_inf = zeros(1, n);
qlim_sup = zeros(1, n);
for j = 1:n
    qlim_inf(j) = robot.links(j).qlim(1);
    qlim_sup(j) = robot.links(j).qlim(2);
end

for intento = 1:N_intentos
    q_rand = qlim_inf + rand(1, n) .* (qlim_sup - qlim_inf);
    q_sol = robot.ikine(T, 'q0', q_rand, 'mask', mascara);

    if isempty(q_sol) || any(isnan(q_sol))
        continue
    end

    if ~dentro_de_limites(robot, q_sol)
        continue
    end
    %Verificamos si existe una solucion similar
    es_nueva = true;
    for k = 1:size(soluciones, 1)
        if norm(q_sol - soluciones(k,:)) < tol
            es_nueva = false;
            break
        end
    end

    if es_nueva
        soluciones(end+1, :) = q_sol; %#ok<AGROW>
    end
end
end


function resultado_limpio = eliminar_repetidas(ramas, tol)
%Eliminamos trayectorias completas que sean casi identicas
%Usamos la norma de frobenius para comparar matrices de trayectoria
n_ramas   = length(ramas);
conservar = true(n_ramas, 1);

for i = 1:n_ramas
    if ~conservar(i), continue; end
    for j = i+1:n_ramas
        if ~conservar(j), continue; end
        if norm(ramas{i} - ramas{j}, 'fro') < tol * size(ramas{i}, 1)
            conservar(j) = false;
        end
    end
end

resultado_limpio = ramas(conservar);
end