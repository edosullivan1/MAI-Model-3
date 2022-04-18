close all;clear;
%USING NEWTON-EULER MECHANICS
%https://www.mdpi.com/1424-8220/15/5/11258/htm#b17-sensors-15-11258

%24-27 seconds deep squat
pos = csvread('VICON_MASTER.csv',4808,0,'A4809..AII5409');

%GRF = csvread('VICON_MASTER.csv',38863,0,'A38864..N41864');
%GRFz_L = GRF(:,4);
%GRFz_R = GRF(:,10);
%10-13 seconds deep squat
%pos = csvread('VICON_MASTER.csv',2008,0,'A2009..AII2609');
Vicon_Moment = readmatrix('VICON_MOMENT_MAG_KG_ADJUSTED.xlsx','Range','A4:AG604');


% Remember that these parameters are zero-based, 
% so that column A maps to 0 and row 252 maps to 251.
t = pos(:, 1);

Ax = pos(:, 101); Ay = pos(:, 102); Az = pos(:, 103);   %Left Toe
Bx = pos(:, 380); By = pos(:, 381); Bz = pos(:, 382);   %Left Ankle
Cx = pos(:, 368); Cy = pos(:, 369); Cz = pos(:, 370);   %Left Knee
Dx = pos(:, 377); Dy = pos(:, 378); Dz = pos(:, 379);   %Left Hip
Ex = pos(:, 521); Ey = pos(:, 522); Ez = pos(:, 523);   %Left Shoulder

Fx = pos(:, 14); Fy = pos(:, 15); Fz = pos(:, 16);    %Neck
Gx = pos(:, 464); Gy = pos(:, 465); Gz = pos(:, 466);    %Head

Hx = pos(:, 125); Hy = pos(:, 126); Hz = pos(:, 127);   %Right Toe
Ix = pos(:, 428); Iy = pos(:, 429); Iz = pos(:, 430);   %Right Ankle
Jx = pos(:, 416); Jy = pos(:, 417); Jz = pos(:, 418);   %Right Knee
Kx = pos(:, 425); Ky = pos(:, 426); Kz = pos(:, 427);   %Right Hip
Lx = pos(:, 557); Ly = pos(:, 558); Lz = pos(:, 559);   %Right Shoulder
Mx = pos(:, 356); My = pos(:, 357); Mz = pos(:, 358);   %Pelvis
Nx = pos(:, 524); Ny = pos(:, 525); Nz = pos(:, 526);   %Left Wrist
Ox = pos(:, 560); Oy = pos(:, 561); Oz = pos(:, 562);   %Right Wrist
Px = pos(:, 98); Py = pos(:, 99); Pz = pos(:, 100);   %Left Heel
Qx = pos(:, 122); Qy = pos(:, 123); Qz = pos(:, 124);   %Right Heel



Sx = [Ax, Bx, Cx, Dx, Ex, Nx, Ox, Fx, Lx, Kx, Mx, Gx, Hx, Ix, Jx];
Sy = [Ay, By, Cy, Dy, Ey, Ny, Mz, Fy, Ly, Ky, My, Gy, Hy, Iy, Jy];
Sz = [Az, Bz, Cz, Dz, Ez, Nz, Oz, Fz, Lz, Kz, Mz, Gz, Hz, Iz, Jz];


Vicon_t = Vicon_Moment(:,1);

M_Vicon_Hip_L = Vicon_Moment(:,17);
M_Vicon_Hip_R = Vicon_Moment(:,22);

M_Vicon_Knee_L = Vicon_Moment(:,6);
M_Vicon_Knee_R = Vicon_Moment(:,11);

M_Vicon_Ankle_L = Vicon_Moment(:,28);
M_Vicon_Ankle_R = Vicon_Moment(:,33);

m_total = 80;
L_total = 1.86;

density = 1000;
g = [0, 0, 9.81];
g_x = 0;
g_y = 0;
g_z = -9.81;

head_scale = 1.1;

% Estimates for trunk are from Plagenhoef et al., 1983 
% https://exrx.net/Kinesiology/Segments

m_foot_L  = 0.5*0.0143*m_total;
m_leg_L   = 0.5*0.0475*m_total;
m_thigh_L = 0.5*0.105*m_total;
m_pelvis_L = 0.5*0.1366*m_total;  %w/0 pelvis
m_arm_L = 0.5*0.057*m_total;

m_torso = 0.551*m_total;     %not including arms

m_foot_R  = 0.5*0.0143*m_total;
m_leg_R   = 0.5*0.0475*m_total;
m_thigh_R = 0.5*0.105*m_total;    %w/o pelvis
m_pelvis_R = 0.5*0.1366*m_total;  

m_arm_R = 0.5*0.057*m_total;

m_pelvis = 0.1366*m_total; 
m_head_neck  = 0.0826*m_total;

m_torso_arms = m_arm_R + m_arm_L + m_torso;


num_rows = height(Ax);

for i = 1:num_rows

L_foot_L(i) = sqrt( (Ax(i) - Bx(i))^2 + (Ay(i) - By(i))^2 + (Az(i) - Bz(i))^2);
L_leg_L(i) = sqrt( (Cx(i) - Bx(i))^2 + (Cy(i) - By(i))^2 + (Cz(i) - Bz(i))^2);
L_thigh_L(i) = sqrt( (Dx(i) - Cx(i))^2 + (Dy(i) - Cy(i))^2 + (Dz(i) - Cz(i))^2);
L_torso_L(i) = sqrt( (Fx(i) - Mx(i))^2 + (Fy(i) - My(i))^2 + (Fz(i) - Mz(i))^2);
L_arm_L(i) = sqrt( (Nx(i) - Ex(i))^2 + (Ny(i) - Ey(i))^2 + (Nz(i) - Ez(i))^2);

L_foot_R(i) = sqrt( (Hx(i) - Ix(i))^2 + (Hy(i) - Iy(i))^2 + (Hz(i) - Iz(i))^2);
L_leg_R(i) = sqrt( (Jx(i) - Ix(i))^2 + (Jy(i) - Iy(i))^2 + (Jz(i) - Iz(i))^2);
L_thigh_R(i) = sqrt( (Kx(i) - Jx(i))^2 + (Ky(i) - Jy(i))^2 + (Kz(i) - Jz(i))^2);
L_torso_R(i) = sqrt( (Fx(i) - Mx(i))^2 + (Fy(i) - My(i))^2 + (Fz(i) - Mz(i))^2);
L_arm_R(i) = sqrt( (Ox(i) - Lx(i))^2 + (Oy(i) - Ly(i))^2 + (Oz(i) - Lz(i))^2);

L_head_neck(i)  = head_scale*sqrt( (Gx(i) - Fx(i))^2 + (Gy(i) - Fy(i))^2 + (Gz(i) - Fz(i))^2);

L_check(i) = L_leg_L(i) + L_thigh_L(i) + L_torso_L(i) + L_head_neck(i);

theta_foot_L(i)  = atand((By(i) - Ay(i))/sqrt(((Bx(i) - Ax(i))^2) + ((Bz(i) - Az(i))^2)));
theta_leg_L(i)   = atand((Cy(i) - By(i))/sqrt(((Cx(i) - Bx(i))^2) + ((Cz(i) - Bz(i))^2)));
theta_thigh_L(i) = atand((Dy(i) - Cy(i))/sqrt(((Dx(i) - Cx(i))^2) + ((Dz(i) - Cz(i))^2)));
theta_arm_L(i) = atand((Ey(i) - Ny(i))/sqrt(((Ex(i) - Nx(i))^2) + ((Ez(i) - Nz(i))^2)));

theta_torso(i) = atand((Fy(i) - My(i))/sqrt(((Fx(i) - Mx(i))^2) + ((Fz(i) - Mz(i))^2)));

theta_foot_R(i)  = atand((Iy(i) - Hy(i))/sqrt(((Ix(i) - Hx(i))^2) + ((Iz(i) - Hz(i))^2)));
theta_leg_R(i) = atand((Jy(i) - Iy(i))/sqrt(((Jx(i) - Ix(i))^2) + ((Jz(i) - Iz(i))^2)));
theta_thigh_R(i) = atand((Ky(i) - Jy(i))/sqrt(((Kx(i) - Jx(i))^2) + ((Kz(i) - Jz(i))^2)));
theta_arm_R(i) = atand((Ly(i) - Oy(i))/sqrt(((Lx(i) - Ox(i))^2) + ((Lz(i) - Oz(i))^2)));

theta_head(i)  = atand((Gy(i) - Fy(i))/sqrt(((Gx(i) - Fx(i))^2) + ((Gz(i) - Fz(i))^2)));


foot_L_cg_x(i) = (Ax(i) + Bx(i))/2;
foot_L_cg_y(i) = (Ay(i) + By(i))/2;
foot_L_cg_z(i) = (Az(i) + Bz(i))/2;

leg_L_cg_x(i) =  (Bx(i) + Cx(i))/2;
leg_L_cg_y(i) =  (By(i) + Cy(i))/2;
leg_L_cg_z(i) =  (Bz(i) + Cz(i))/2;

thigh_L_cg_x(i) =  (Cx(i) + Dx(i))/2;
thigh_L_cg_y(i) =  (Cy(i) + Dy(i))/2;
thigh_L_cg_z(i) =  (Cz(i) + Dz(i))/2;

arm_L_cg_x(i) =  (Nx(i) + Ex(i))/2;
arm_L_cg_y(i) =  (Ny(i) + Ey(i))/2;
arm_L_cg_z(i) =  (Nz(i) + Ez(i))/2;

torso_L_cg_x(i) =  (Ex(i) + Dx(i))/2;
torso_L_cg_y(i) =  (Ey(i) + Dy(i))/2;
torso_L_cg_z(i) =  (Ez(i) + Dz(i))/2;

pelvis_L_cg_x(i) = (Mx(i) + Dx(i))/2;
pelvis_L_cg_y(i) =  (My(i) + Dy(i))/2;
pelvis_L_cg_z(i) =  (Mz(i) + Dz(i))/2;

torso_cg_x(i) =  (Fx(i) + Mx(i))/2;
torso_cg_y(i) =  (Fy(i) + My(i))/2;     %Middle torso
torso_cg_z(i) =  (Fz(i) + Mz(i))/2;

torso_R_cg_x(i) =  (Lx(i) + Kx(i))/2;
torso_R_cg_y(i) =  (Ly(i) + Ky(i))/2;
torso_R_cg_z(i) =  (Lz(i) + Kz(i))/2;

head_neck_cg_x(i) =  (Fx(i) + Gx(i))/2;
head_neck_cg_y(i) =  (Fy(i) + Gy(i))/2;
head_neck_cg_z(i) =  (head_scale * (Fz(i) + Gz(i)))/2;

foot_R_cg_x(i) = (Ix(i) + Hx(i))/2;
foot_R_cg_y(i) = (Iy(i) + Hy(i))/2;
foot_R_cg_z(i) = (Iz(i) + Hz(i))/2;

leg_R_cg_x(i) =  (Jx(i) + Ix(i))/2;
leg_R_cg_y(i) =  (Jy(i) + Iy(i))/2;
leg_R_cg_z(i) =  (Jz(i) + Iz(i))/2;

thigh_R_cg_x(i) =  (Kx(i) + Jx(i))/2;
thigh_R_cg_y(i) =  (Ky(i) + Jy(i))/2;
thigh_R_cg_z(i) =  (Kz(i) + Jz(i))/2;

arm_R_cg_x(i) =  (Ox(i) + Lx(i))/2;
arm_R_cg_y(i) =  (Oy(i) + Ly(i))/2;
arm_R_cg_z(i) =  (Oz(i) + Lz(i))/2;

pelvis_R_cg_x(i) = (Mx(i) + Kx(i))/2;
pelvis_R_cg_y(i) =  (My(i) + Ky(i))/2;
pelvis_R_cg_z(i) =  (Mz(i) + Kz(i))/2;


% Whole body cg


%F_z_ankle_R(i) = GRFz_R(i + 4*(i - 1)) - m_foot_R * g(3);
%F_z_ankle_L(i) = GRFz_L(i + 4*(i - 1)) - m_foot_L * g(3);

%Vector r direction (cg - joint centre)
%GRF acting through cg of foot
%r_1_L(i,:) = [(foot_L_cg_x(i) - ), (foot_L_cg_y(i) - Py(i)), (foot_R_cg_z(i) - Pz(i))];
%r_1_R(i,:) = [(foot_R_cg_x(i) - Qx(i)), (foot_R_cg_y(i) - Qy(i)), (foot_R_cg_z(i) - Qz(i))];

r_2_L(i,:) = [(foot_L_cg_x(i) - Bx(i)), (foot_L_cg_y(i) - By(i)), (foot_L_cg_z(i) - Bz(i))];
r_2_R(i,:) = [(foot_R_cg_x(i) - Ix(i)), (foot_R_cg_y(i) - Iy(i)), (foot_R_cg_z(i) - Iz(i))];

r_3_L(i,:) = [(leg_L_cg_x(i) - Bx(i)), (leg_L_cg_y(i) - By(i)), (leg_L_cg_z(i) - Bz(i))];
r_3_R(i,:) = [(leg_R_cg_x(i) - Ix(i)), (leg_R_cg_y(i) - Iy(i)), (leg_R_cg_z(i) - Iz(i))];

r_4_L(i,:) = [(leg_L_cg_x(i) - Cx(i)), (leg_L_cg_y(i) - Cy(i)), (leg_L_cg_z(i) - Cz(i))];
r_4_R(i,:) = [(leg_R_cg_x(i) - Jx(i)), (leg_R_cg_y(i) - Jy(i)), (leg_R_cg_z(i) - Jz(i))];

r_5_L(i,:) = [(thigh_L_cg_x(i) - Cx(i)), (thigh_L_cg_y(i) - Cy(i)), (thigh_L_cg_z(i) - Cz(i))];
r_5_R(i,:) = [(thigh_R_cg_x(i) - Jx(i)), (thigh_R_cg_y(i) - Jy(i)), (thigh_R_cg_z(i) - Jz(i))];

r_6_L(i,:) = [(thigh_L_cg_x(i) - Dx(i)), (thigh_L_cg_y(i) - Dy(i)), (thigh_L_cg_z(i) - Dz(i))];
r_6_R(i,:) = [(thigh_R_cg_x(i) - Kx(i)), (thigh_R_cg_y(i) - Ky(i)), (thigh_R_cg_z(i) - Kz(i))];

r_7_L(i,:) = [(Mx(i) - Dx(i)), (My(i) - Dy(i)), (Mz(i) - Dz(i))];
r_7_R(i,:) = [(Mx(i) - Kx(i)), (My(i) - Ky(i)), (Mz(i) - Kz(i))];

r_8(i,:) = [(torso_cg_x(i) - Mx(i)), (torso_cg_y(i) - My(i)), (torso_cg_z(i) - Mz(i))];
r_9(i,:) = [(torso_cg_x(i) - Fx(i)), (torso_cg_y(i) - Fy(i)), (torso_cg_z(i) - Fz(i))];

r_10(i,:) = [(head_neck_cg_x(i)- Fx(i)),(head_neck_cg_y(i)-Fy(i)), (head_neck_cg_z(i)-Fz(i))];

I_ankle_L = 0;
w_dot_L = 0;
I_leg_L = 0;
w_dot_leg_L = 0;

%A(1,1) = 1;
%A(1,3) = -1;
%A(2,2) = 1;
%A(2,4) = -1;
%A(3,3) = 1;
%A(3,5) = -1;
%A(4,4) = 1;
%A(4,6) = -1;
%A(5,5) = 1;
%A(5,7) = -1;
%A(6,6) = 1;
%A(6,8) = -1;
%A(7,7) = 1;
%A(8,8) = 1;

A = [1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	r_2_L(i,3)	-r_2_L(i,2)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	-r_2_L(i,3)	0	r_2_L(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	r_2_L(i,2)	-r_2_L(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	r_2_R(i,3)	-r_2_R(i,2)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	-r_2_R(i,3)	0	r_2_R(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	r_2_R(i,2)	-r_2_R(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	-r_3_L(i,3)	r_3_L(i,2)	0	0	0	0	r_4_L(i,3)	-r_4_L(i,2)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	r_3_L(i,3)	0	-r_3_L(i,1)	0	0	0	-r_4_L(i,3)	0	r_4_L(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	-r_3_L(i,2)	r_3_L(i,1)	0	0	0	0	r_4_L(i,2)	-r_4_L(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	-r_3_R(i,3)	r_3_R(i,2)	0	0	0	0	r_4_R(i,3)	-r_4_R(i,2)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	r_3_R(i,3)	0	-r_3_R(i,1)	0	0	0	-r_4_R(i,3)	0	r_4_R(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	-r_3_R(i,2)	r_3_R(i,1)	0	0	0	0	r_4_R(i,2)	-r_4_R(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	-r_5_L(i,3)	r_5_L(i,2)	0	0	0	0	r_6_L(i,3)	-r_6_L(i,2)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	r_5_L(i,3)	0	-r_5_L(i,1)	0	0	0	-r_6_L(i,3)	0	r_6_L(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	-r_5_L(i,2)	r_5_L(i,1)	0	0	0	0	r_6_L(i,2)	-r_6_L(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-r_5_R(i,3)	r_5_R(i,2)	0	0	0	0	r_6_R(i,3)	-r_6_R(i,2)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	r_5_R(i,3)	0	-r_5_R(i,1)	0	0	0	-r_6_R(i,3)	0	r_6_R(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-r_5_R(i,2)	r_5_R(i,1)	0	0	0	0	r_6_R(i,2)	-r_6_R(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	-1	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-r_7_L(i,3)	r_7_L(i,2)	0	r_7_R(i,3)	-r_7_R(i,2)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	-1	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	r_7_L(i,3)	0	-r_7_L(i,1)	-r_7_R(i,3)	0	r_7_R(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	-1	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-r_7_L(i,2)	r_7_L(i,1)	0	r_7_R(i,2)	-r_7_R(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	-1	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-r_8(i,3)	r_8(i,2)	0	r_9(i,3)	-r_9(i,2)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	-1	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	r_8(i,3)	0	-r_8(i,1)	-r_9(i,3)	0	r_9(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	-1	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-r_8(i,2)	r_8(i,1)	0	r_9(i,2)	-r_9(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	-1	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-r_10(i,3)	r_10(i,2)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	r_10(i,3)	0	-r_10(i,1)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-r_10(i,2)	r_10(i,2)	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	;
];

B = zeros(54,1);
B(3) = m_foot_L * g_z; 
B(6) = m_foot_R * g_z; 
B(9) = m_leg_L * g_z; 
B(12) = m_leg_R * g_z;
B(15) = m_thigh_L * g_z; 
B(18) = m_thigh_R * g_z; 
B(21) = m_pelvis * g_z; 
B(24) = m_torso_arms * g_z;
B(51) = m_head_neck * g_z;

X = pinv(A)*B;
%X = A\B;


F_GRF_L_x(i) = X(1);
F_GRF_L_y(i) = X(2);
F_GRF_L_z(i) = X(3);
F_GRF_R_x(i) = X(4);
F_GRF_R_y(i) = X(5);
F_GRF_R_z(i) = X(6);
F_ANKLE_L_x(i) = X(7);
F_ANKLE_L_y(i) = X(8);
F_ANKLE_L_z(i) = X(9);
F_ANKLE_R_x(i) = X(10);
F_ANKLE_R_y(i) = X(11);
F_ANKLE_R_z(i) = X(12);
F_KNEE_L_x(i) = X(13);
F_KNEE_L_y(i) = X(14);
F_KNEE_L_z(i) = X(15);
F_KNEE_R_x(i) = X(16);
F_KNEE_R_y(i) = X(17);
F_KNEE_R_z(i) = X(18);
F_HIP_L_x(i) = X(19);
F_HIP_L_y(i) = X(20);
F_HIP_L_z(i) = X(21);
F_HIP_R_x(i) = X(22);
F_HIP_R_y(i) = X(23);
F_HIP_R_z(i) = X(24);
F_PELVIS_x(i) = X(25);
F_PELVIS_y(i) = X(26);
F_PELVIS_z(i) = X(27);
F_HEAD_NECK_x(i) = X(28);
F_HEAD_NECK_y(i) = X(29);
F_HEAD_NECK_z(i) = X(30);
M_Ankle_L_x(i) = X(31);
M_Ankle_L_y(i) = X(32);
M_Ankle_L_z(i) = X(33);
M_Ankle_R_x(i) = X(34);
M_Ankle_R_y(i) = X(35);
M_Ankle_R_z(i) = X(36);
M_Knee_L_x(i) = X(37);
M_Knee_L_y(i) = X(38);
M_Knee_L_z(i) = X(39);
M_Knee_R_x(i) = X(40);
M_Knee_R_y(i) = X(41);
M_Knee_R_z(i) = X(42);
M_Hip_L_x(i) = X(43);
M_Hip_L_y(i) = X(44);
M_Hip_L_z(i) = X(45);
M_Hip_R_x(i) = X(46);
M_Hip_R_y(i) = X(47);
M_Hip_R_z(i) = X(48);
M_Pelvis_x(i) = X(49);
M_Pelvis_y(i) = X(50);
M_Pelvis_z(i) = X(51);
M_Head_Neck_x(i) = X(52);
M_Head_Neck_y(i) = X(53);
M_Head_Neck_z(i) = X(54);



M_Hip_L_mag(i) = sqrt((M_Hip_L_x(i).^2) + (M_Hip_L_y(i).^2) + (M_Hip_L_z(i).^2));
M_Hip_R_mag(i) = sqrt((M_Hip_R_x(i).^2) + (M_Hip_R_y(i).^2) + (M_Hip_R_z(i).^2));

M_Knee_L_mag(i) = sqrt((M_Knee_L_x(i).^2) + (M_Knee_L_y(i).^2) + (M_Knee_L_z(i).^2));
M_Knee_R_mag(i) = sqrt((M_Knee_R_x(i).^2) + (M_Knee_R_y(i).^2) + (M_Knee_R_z(i).^2));

M_Ankle_L_mag(i) = sqrt((M_Ankle_L_x(i).^2) + (M_Ankle_L_y(i).^2) + (M_Ankle_L_z(i).^2));
M_Ankle_R_mag(i) = sqrt((M_Ankle_R_x(i).^2) + (M_Ankle_R_y(i).^2) + (M_Ankle_R_z(i).^2));

%Check
%Equation 1
Eq_1_R_x = F_GRF_R_x(i) - F_ANKLE_R_x(i);

%Equation 6 L
Eq_6_L_x(i) = - M_Ankle_L_x(i) - r_2_L(i,2)*F_ANKLE_L_z(i) + r_2_L(i,3)*F_ANKLE_L_y(i); 

Eq_6_L_y(i) = - M_Ankle_L_y(i) - r_2_L(i,3)*F_ANKLE_L_x(i) + r_2_L(i,1)*F_ANKLE_L_z(i);

Eq_6_L_z(i) = - M_Ankle_L_z(i) - r_2_L(i,1)*F_ANKLE_L_y(i) + r_2_L(i,2)*F_ANKLE_L_x(i);

%Equation 6 R
Eq_6_R_x(i) = - M_Ankle_R_x(i) - r_2_R(i,2)*F_ANKLE_R_z(i) + r_2_R(i,3)*F_ANKLE_R_y(i); 

Eq_6_R_y(i) = - M_Ankle_R_y(i) - r_2_R(i,3)*F_ANKLE_R_x(i) + r_2_R(i,1)*F_ANKLE_R_z(i);

Eq_6_R_z(i) = - M_Ankle_R_z(i) - r_2_R(i,1)*F_ANKLE_R_y(i) + r_2_R(i,2)*F_ANKLE_R_x(i);

%Equation 7 L
Eq_7_L_x(i) = M_Ankle_L_x(i) -M_Knee_L_x(i) ...
    + r_3_L(i,2)*F_ANKLE_L_z(i) - r_3_L(i,3)*F_ANKLE_L_y(i)...
    - r_4_L(i,2)*F_KNEE_L_z(i) + r_4_L(i,3)*F_KNEE_L_y(i);

Eq_7_L_y(i) = M_Ankle_L_y(i) -M_Knee_L_y(i) ...
    + r_3_L(i,3)*F_ANKLE_L_x(i) - r_3_L(i,1)*F_ANKLE_L_z(i)...
    - r_4_L(i,3)*F_KNEE_L_x(i) + r_4_L(i,1)*F_KNEE_L_z(i);

Eq_7_L_z(i) = M_Ankle_L_z(i) -M_Knee_L_z(i) ...
    + r_3_L(i,1)*F_ANKLE_L_y(i) - r_3_L(i,2)*F_ANKLE_L_x(i)...
    - r_4_L(i,1)*F_KNEE_L_y(i) + r_4_L(i,2)*F_KNEE_L_x(i);

%Equation 7 R
Eq_7_R_x(i) = M_Ankle_R_x(i) -M_Knee_R_x(i) ...
    + r_3_R(i,2)*F_ANKLE_R_z(i) - r_3_R(i,3)*F_ANKLE_R_y(i)...
    - r_4_R(i,2)*F_KNEE_R_z(i) + r_4_R(i,3)*F_KNEE_R_y(i);

Eq_7_R_y(i) = M_Ankle_R_y(i) -M_Knee_R_y(i) ...
    + r_3_R(i,3)*F_ANKLE_R_x(i) - r_3_R(i,1)*F_ANKLE_R_z(i)...
    - r_4_R(i,3)*F_KNEE_R_x(i) + r_4_R(i,1)*F_KNEE_R_z(i);

Eq_7_R_z(i) = M_Ankle_R_z(i) -M_Knee_R_z(i) ...
    + r_3_R(i,1)*F_ANKLE_R_y(i) - r_3_R(i,2)*F_ANKLE_R_x(i)...
    - r_4_R(i,1)*F_KNEE_R_y(i) + r_4_R(i,2)*F_KNEE_R_x(i);
end


figure (1)
plot(t,M_Knee_L_mag)
title('Torque-Time plot of Left Knee')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t,M_Vicon_Knee_L)

figure (2)
plot(t,M_Knee_R_mag)
title('Torque-Time plot of Right Knee')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t,M_Vicon_Knee_R)

figure (3)
plot(t,M_Hip_L_mag)
title('Torque-Time plot of Left Hip')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t,M_Vicon_Hip_L)

figure (4)
plot(t,M_Hip_R_mag)
title('Torque-Time plot of Right Hip')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t,M_Vicon_Hip_R)

figure (5)
plot(t,M_Ankle_L_mag)
title('Torque-Time plot of Left Ankle')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t,M_Vicon_Ankle_L)


figure (6)
plot(t,M_Ankle_R_mag)
title('Torque-Time plot of Right Ankle')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t,M_Vicon_Ankle_R)

CGx = [foot_L_cg_x;leg_L_cg_x;thigh_L_cg_x;torso_cg_x;head_neck_cg_x;thigh_R_cg_x;leg_R_cg_x;foot_R_cg_x];
CGy = [foot_L_cg_y;leg_L_cg_y;thigh_L_cg_y;torso_cg_y;head_neck_cg_y;thigh_R_cg_y;leg_R_cg_y;foot_R_cg_y];
CGz = [foot_L_cg_z;leg_L_cg_z;thigh_L_cg_z;torso_cg_z;head_neck_cg_z;thigh_R_cg_z;leg_R_cg_z;foot_R_cg_z];


figure(7)
plot3(Sx(1,:), Sy(1,:), Sz(1,:), 'o', 'Color', 'b', 'MarkerSize',4,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,1), CGy(:,1), CGz(:,1), 'o', 'Color', 'g', 'MarkerSize',4,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
hold on
quiver3(Bx(1), By(1), Bz(1), r_2_L(1,1), r_2_L(1,2), r_2_L(1,3), 'LineWidth',3, 'Color','g', 'Marker', '*', 'MarkerSize',6)
quiver3(Bx(1), By(1), Bz(1), r_3_L(1,1), r_3_L(1,2), r_3_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Cx(1), Cy(1), Cz(1), r_4_L(1,1), r_4_L(1,2), r_4_L(1,3), 'LineWidth',3, 'Color','c', 'Marker', '*', 'MarkerSize',6)
quiver3(Cx(1), Cy(1), Cz(1), r_5_L(1,1), r_5_L(1,2), r_5_L(1,3), 'LineWidth',3, 'Color','m', 'Marker', '*', 'MarkerSize',6)
quiver3(Dx(1), Dy(1), Dz(1), r_6_L(1,1), r_6_L(1,2), r_6_L(1,3), 'LineWidth',3, 'Color','y', 'Marker', '*', 'MarkerSize',6)
quiver3(Dx(1), Dy(1), Dz(1), r_7_L(1,1), r_7_L(1,2), r_7_L(1,3), 'LineWidth',3, 'Color','k', 'Marker', '*', 'MarkerSize',6)
quiver3(Ix(1), Iy(1), Iz(1), r_2_R(1,1), r_2_R(1,2), r_2_R(1,3), 'LineWidth',3, 'Color','g', 'Marker', '*', 'MarkerSize',6)
quiver3(Ix(1), Iy(1), Iz(1), r_3_R(1,1), r_3_R(1,2), r_3_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Jx(1), Jy(1), Jz(1), r_4_R(1,1), r_4_R(1,2), r_4_R(1,3), 'LineWidth',3, 'Color','c', 'Marker', '*', 'MarkerSize',6)
quiver3(Jx(1), Jy(1), Jz(1), r_5_R(1,1), r_5_R(1,2), r_5_R(1,3), 'LineWidth',3, 'Color','m', 'Marker', '*', 'MarkerSize',6)
quiver3(Kx(1), Ky(1), Kz(1), r_6_R(1,1), r_6_R(1,2), r_6_R(1,3), 'LineWidth',3, 'Color','y', 'Marker', '*', 'MarkerSize',6)
quiver3(Kx(1), Ky(1), Kz(1), r_7_R(1,1), r_7_R(1,2), r_7_R(1,3), 'LineWidth',3, 'Color','k', 'Marker', '*', 'MarkerSize',6)
quiver3(Mx(1), My(1), Mz(1), r_8(1,1), r_8(1,2), r_8(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Fx(1), Fy(1), Fz(1), r_9(1,1), r_9(1,2), r_9(1,3), 'LineWidth',3, 'Color','g', 'Marker', '*', 'MarkerSize',6)
quiver3(Fx(1), Fy(1), Fz(1), r_10(1,1), r_10(1,2), r_10(1,3), 'LineWidth',3, 'Color','r', 'Marker', '*', 'MarkerSize',6)

axis equal

% figure(8)
% plot3(Sx(200,:), Sy(200,:), Sz(200,:), 'o', 'Color', 'b', 'MarkerSize',4,'MarkerFaceColor','#D9FFFF')
% hold on
% plot3(CGx(:,200), CGy(:,200), CGz(:,200), 'o', 'Color', 'g', 'MarkerSize',4,'MarkerFaceColor','#D9FFFF')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% quiver3(Bx(200), By(200), Bz(200), r_2_L(200,1), r_2_L(200,2), r_2_L(200,3), 'LineWidth',3, 'Color','g', 'Marker', '*', 'MarkerSize',6)
% quiver3(Bx(200), By(200), Bz(200), r_3_L(200,1), r_3_L(200,2), r_3_L(200,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
% quiver3(Cx(200), Cy(200), Cz(200), r_4_L(200,1), r_4_L(200,2), r_4_L(200,3), 'LineWidth',3, 'Color','c', 'Marker', '*', 'MarkerSize',6)
% quiver3(Cx(200), Cy(200), Cz(200), r_5_L(200,1), r_5_L(200,2), r_5_L(200,3), 'LineWidth',3, 'Color','m', 'Marker', '*', 'MarkerSize',6)
% quiver3(Dx(200), Dy(200), Dz(200), r_6_L(200,1), r_6_L(200,2), r_6_L(200,3), 'LineWidth',3, 'Color','y', 'Marker', '*', 'MarkerSize',6)
% quiver3(Dx(200), Dy(200), Dz(200), r_7_L(200,1), r_7_L(200,2), r_7_L(200,3), 'LineWidth',3, 'Color','k', 'Marker', '*', 'MarkerSize',6)
% quiver3(Ix(200), Iy(200), Iz(200), r_2_R(200,1), r_2_R(200,2), r_2_R(200,3), 'LineWidth',3, 'Color','g', 'Marker', '*', 'MarkerSize',6)
% quiver3(Ix(200), Iy(200), Iz(200), r_3_R(200,1), r_3_R(200,2), r_3_R(200,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
% quiver3(Jx(200), Jy(200), Jz(200), r_4_R(200,1), r_4_R(200,2), r_4_R(200,3), 'LineWidth',3, 'Color','c', 'Marker', '*', 'MarkerSize',6)
% quiver3(Jx(200), Jy(200), Jz(200), r_5_R(200,1), r_5_R(200,2), r_5_R(200,3), 'LineWidth',3, 'Color','m', 'Marker', '*', 'MarkerSize',6)
% quiver3(Kx(200), Ky(200), Kz(200), r_6_R(200,1), r_6_R(200,2), r_6_R(200,3), 'LineWidth',3, 'Color','y', 'Marker', '*', 'MarkerSize',6)
% quiver3(Kx(200), Ky(200), Kz(200), r_7_R(200,1), r_7_R(200,2), r_7_R(200,3), 'LineWidth',3, 'Color','k', 'Marker', '*', 'MarkerSize',6)
% quiver3(Mx(200), My(200), Mz(200), r_8(200,1), r_8(200,2), r_8(200,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
% quiver3(Fx(200), Fy(200), Fz(200), r_9(200,1), r_9(200,2), r_9(200,3), 'LineWidth',3, 'Color','g', 'Marker', '*', 'MarkerSize',6)
% quiver3(Fx(200), Fy(200), Fz(200), r_10(200,1), r_10(200,2), r_10(200,3), 'LineWidth',3, 'Color','r', 'Marker', '*', 'MarkerSize',6)
% axis equal

figure(9)
plot(Eq_6_L_x)
title('Equation 6 L x')
xlabel('i')
ylabel('Deviation')

figure(10)
plot(Eq_6_L_y)
title('Equation 6 L y')
xlabel('i')
ylabel('Deviation')

figure(11)
plot(Eq_6_L_z)
title('Equation 6 L z')
xlabel('i')
ylabel('Deviation')

figure(12)
plot(Eq_6_R_x)
title('Equation 6 R x')
xlabel('i')
ylabel('Deviation')

figure(13)
plot(Eq_6_R_y)
title('Equation 6 R y')
xlabel('i')
ylabel('Deviation')

figure(14)
plot(Eq_6_R_z)
title('Equation 6 R z')
xlabel('i')
ylabel('Deviation')

figure(15)
plot(Eq_7_L_x)
title('Equation 7 L x')
xlabel('i')
ylabel('Deviation')

figure(16)
plot(Eq_7_L_y)
title('Equation 7 L y')
xlabel('i')
ylabel('Deviation')

figure(17)
plot(Eq_7_L_z)
title('Equation 7 L z')
xlabel('i')
ylabel('Deviation')

figure(18)
plot(Eq_7_R_x)
title('Equation 7 R x')
xlabel('i')
ylabel('Deviation')

figure(19)
plot(Eq_7_R_y)
title('Equation 7 R y')
xlabel('i')
ylabel('Deviation')

figure(20)
plot(Eq_7_R_z)
title('Equation 7 R z')
xlabel('i')
ylabel('Deviation')


figure (21)
subplot(3, 2, 1)
plot(t,M_Hip_L_mag)
title('Torque-Time plot of Left Hip')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Hip_L)
ylim([0 8e4])

subplot(3,2,2)
plot(t,M_Hip_R_mag)
title('Torque-Time plot of Right Hip')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Hip_R)
ylim([0 8e4])

subplot(3,2,3)
plot(t,M_Knee_L_mag)
title('Torque-Time plot of Left Knee')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Knee_L)
ylim([0 12e4])

subplot(3,2,4)
plot(t,M_Knee_R_mag)
title('Torque-Time plot of Right Knee')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Knee_R)
ylim([0 12e4])

subplot(3,2,5)
plot(t,M_Ankle_L_mag)
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Ankle_L)
title('Torque-Time plot of Left Ankle')
ylim([0 4e4])

subplot(3,2,6)
plot(t, M_Ankle_R_mag)
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Ankle_R)
title('Torque-Time plot of Right Ankle')
ylim([0 4e4])

figure(22)
subplot(3,2,1)
plot(Eq_6_L_x)
title('Equation 3.7 Left Foot (x-component)')
xlabel('i')
ylabel('Deviation (Nmm)')
xlim([0 601])

subplot(3,2,2)
plot(Eq_7_L_x)
title('Equation 3.8 Left Leg (x-component)')
xlabel('i')
ylabel('Deviation (Nmm)')
xlim([0 601])

subplot(3,2,3)
plot(Eq_6_L_y)
title('Equation 3.7 Left Foot (y-component')
xlabel('i')
ylabel('Deviation (Nmm)')
xlim([0 601])

subplot(3,2,4)
plot(Eq_7_L_y)
title('Equation 3.8 Left Leg (y-component)')
xlabel('i')
ylabel('Deviation (Nmm)')
xlim([0 601])


subplot(3,2,5)
plot(Eq_6_L_z)
title('Equation 3.7 Left Foot (z-component')
xlabel('i')
ylabel('Deviation (Nmm)')
xlim([0 601])

subplot(3,2,6)
plot(Eq_7_L_z)
title('Equation 3.8 Left Leg (z-component)')
xlabel('i')
ylabel('Deviation (Nmm)')
xlim([0 601])
