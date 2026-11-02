
clear; clc; close all;

% b) Matriz DH con parámetros: tita, d, a, alfa, sigma
% Unidades: metros y radianes 
% Fuente: Datasheet oficial ABB IRB 5720-90/3.0 LID
%Corregimos d4
dh = [
    0, 0.670, 0.388, -pi/2, 0;  % S1 
    0, 0.000, 1.120,     0, 0;  % S2 
    0, 0.000, 0.200, -pi/2, 0;  % S3 
    0, 1.5475, 0.000,  pi/2, 0;  % S4 
    0, 0.000, 0.000, -pi/2, 0;  % S5 
    0, 0.333, 0.000,     0, 0;  % S6 
];

% c) Crear objeto SerialLink
R = SerialLink(dh);

% c.i) Nombre del robot
R.name = 'ABB IRB 5720-90/3.0 LID';

% c.ii) Límites articulares según datasheet (en radianes)
R.links(1).qlim = deg2rad([-170,  170]);  % q1: Eje 1
R.links(2).qlim = deg2rad([ -75,  145]);  % q2: Eje 2
R.links(3).qlim = deg2rad([-160,   70]);  % q3: Eje 3 (LID)
R.links(4).qlim = deg2rad([-300,  300]);  % q4: Eje 4
R.links(5).qlim = deg2rad([-120,  120]);  % q5: Eje 5 (LID)
R.links(6).qlim = deg2rad([-200,  200]);  % q6: Eje 6 (LID)

% c.iii) Offset articular
% El eje 2 en posición cero del fabricante tiene el brazo vertical,
% no horizontal como define DH → se corrige con -pi/2 en q2
R.offset = [0, -pi/2, 0, 0, 0, 0];

% c.iv) Base del robot
% Se eleva 0.3 m del suelo para representar la base
R.base = transl(0, 0, 0.3);

% c.v) Tool: Atomizador rotativo ABB RB1000i
% Longitud total desde brida hasta punta del bell cup: 348 mm = 0.348 m
% El atomizador se monta a lo largo del eje Z del último sistema
R.tool = transl(0, 0, 0.348);

% d) Workspace para ploteo (en metros)
% [-limX, +limX, -limY, +limY, -limZ, +limZ]
workspace = [-4, 4, -4, 4, -0.5, 4];

