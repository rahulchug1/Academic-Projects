
% Assume the bicycle model
% a- loaded case  
W =87800;             % in N
l =3.2;               % wheelbase- in m
a1 =1.9;              % Distance of CG behind front axle -in m
a2 =l-a1;             % Distance of CG infront of back axle-in m
g =9.81;              % in m/s^2
V1 =50;               % in km/h
V =V1*(1000/3600)     % in m/s

% Assume delta_o = delta_i = delta_f 

delta_f =2;           % in degrees (assume)

% Assume the slip angle alfa = 4 degree 

alpha_f =4;
alpha_r = 4;
alpha_f_rad =4*(pi/180);
alpha_r_rad = 4*(pi/180);

% Normal load

Fzf = ((W/2)*(a2/l))
Fzr = ((W/2)*(a1/l))

% from figure 1 and interpolation
% Cornering Forces
Fyf = 8125;        %in N
Fyr = 12785.5;     %in N

C_alpha_f =Fyf/(alpha_f)
C_alpha_r =Fyr/(alpha_r)

C_alpha_f_rad =Fyf/(alpha_f_rad)
C_alpha_r_rad =Fyr/(alpha_r_rad)

K_us_deg =(Fzf/C_alpha_f)-(Fzr/C_alpha_r)
K_us_rad =K_us_deg*(pi/180)

% Yaw Velocity Gain

G_r = V/(l +(K_us_rad*(V*V/g)))
G_r_deg = V/(l +(K_us_deg*(V*V/g)))
% Lateral Acceleration Gain

G_ay=(V^2)/(l*g +(K_us_rad*(V^2)))	
G_ay_deg=(V^2)/(l*g +(K_us_deg*(V^2)))	

% Curvature Gain

G_R = 1/(l+(K_us_rad*(V^2/g)))
G_R_deg = 1/(l+(K_us_deg*(V^2/g)))

% Curvature Radius at staedy state condition

R = (l+K_us_deg*(V^2/g))/delta_f

% curvature radius in stady from transient response
R_t= (l/delta_f)*(1-(((W/9.81)*(a1*C_alpha_f-a2*C_alpha_r)*V^2)/(2*l^2*C_alpha_f*C_alpha_r)))
