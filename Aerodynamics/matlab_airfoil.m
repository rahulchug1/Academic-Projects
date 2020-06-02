
## Input
Vinf = 1;                                                                   % Freestream velocity []
AoA  = 0;                                                                   % Angle of attack [deg]
NACA = '6412';                                                              % NACA airfoil to load [####]

## Convert angle of attack to radians
alpha = AoA*(pi/180);                                                       % Angle of attack [rad]

## Number of boundary points and panels
numPts = length(XB);                                                        % Number of boundary points
numPan = numPts - 1;                                                        % Number of panels (control points)


%% SOURCE PANEL METHOD 

% Initialize variables
XC   = zeros(numPan,1);                                                     % X-coordinate array
YC   = zeros(numPan,1);                                                     % Y-coordinate array
S    = zeros(numPan,1);                                                     % panel length array
phiD = zeros(numPan,1);                                                     % panel orientation angle array [deg]

% geometric quantities for a panel of airfoil
for i = 1:1:numPan                                                          % Loop over all panels
    XC(i)   = 0.5*(XB(i)+XB(i+1));                                          % X-value of control point
    YC(i)   = 0.5*(YB(i)+YB(i+1));                                          % Y-value of control point
    dx      = XB(i+1)-XB(i);                                                % Change in X between boundary points
    dy      = YB(i+1)-YB(i);                                                % Change in Y between boundary points
    S(i)    = (dx^2 + dy^2)^0.5;                                            % Length of the panel
	phiD(i) = atan2d(dy,dx);                                                  % Angle of the panel (positive X-axis to inside face)
    if (phiD(i) < 0)                                                      
        phiD(i) = phiD(i) + 360;
    end
end


deltaD             = phiD + 90;                                             % Angle from positive X-axis to outward normal vector [deg]
betaD              = deltaD - AoA;                                          % Angle between freestream vector and outward normal vector [deg]
betaD(betaD > 360) = betaD(betaD > 360) - 360;                              % Make sure angles aren't greater than 360 [deg]

% Convert angles from [deg] to [rad]
phi  = phiD.*(pi/180);                                                      % in rad
beta = betaD.*(pi/180);                                                     % in rad

%% SOURCE PANEL STRENGTHS CALCULATION

% Geometric integral (normal [I] and tangential [J])

[I,J] = Integral_IJ(XC,YC,XB,YB,phi,S);                                  % Compute geometric integrals

%  Formulation of A matrix

A = zeros(numPan,numPan);                                                   % Initialize the A matrix
for i = 1:1:numPan                                                          % Loop over all i panels
    for j = 1:1:numPan                                                      % Loop over all j panels
        if (i == j)                                                         % If the panels are the same
            A(i,j) = pi;                                                    
        else                                                                % If panels are not the same
            A(i,j) = I(i,j);                                                % Set A equal to geometric integral
        end
    end
end

% FORMULATION OF  b array

b = zeros(numPan,1);                                                        
for i = 1:1:numPan                                                          
    b(i) = -Vinf*2*pi*cos(beta(i));                                         
end

% Compute source panel strengths (lambda) from system of equations
lambda  = A\b;                                                              % source strength values

%  source strenghts sum calculation

sumLambda = sum(lambda.*S);                                                 % Check sum of source panel strengths
fprintf('Sum of L: %g\n',sum(lambda.*S));                                   % Print sum of all source strengths


% Compute velocities

Vt = zeros(numPan,1);                                                       % Initialize tangential velocity array
Cp = zeros(numPan,1);                                                       % Initialize pressure coefficient array
for i = 1:1:numPan                                                          % Loop over all i panels
    addVal  = 0;                                                            % Reset the summation value to zero
    for j = 1:1:numPan                                                      % Loop over all j panels
        addVal = addVal + (lambda(j)/(2*pi))*(J(i,j));                      % Sum all tangential source panel terms
    end
    
    Vt(i) = Vinf*sin(beta(i)) + addVal;                                     % Compute tangential velocity by adding uniform flow term
    Cp(i) = 1-(Vt(i)/Vinf)^2;                                               % Compute pressure coefficient
end

%% COMPUTE CL, CM, CP

M = .06;                                                                      % MP XX 4 digiit airfoil 6412
P = 0.4;
c = 1;                                                                        %chord length
x1 = P;
theta = acos(1-((2*x1)/c))

fun1 = @(th) (M/(P^2))*((2*P)-(2*(c/2)*(1-cos(th))));                         % Integral of DZ/Dx from 0 to theta 
Int_1 = integral(fun1,0,theta);

fun2 = @(th) (M/((1-P)^2))*((2*P)-(2*(c/2)*(1-cos(th))));                     % Integral of DZ/Dx from theta to pi
Int_2 = integral(fun2,theta,pi);
Int = Int_1 + Int_2;                                                          % Integral of DZ/Dx from 0 to pi

A0 = alpha - (Int/pi);                                                        % Ao formula

fun3 = @(th) ((M/(P^2))*((2*P)-(2*(c/2)*(1-cos(th)))).*cos(th));              % Integral of {(DZ/Dx).cos(theta)) w.r.t theta
Int_3 = integral(fun3,0,theta);                                               % limit is 0 to theta 

fun4 = @(th) (M/((1-P)^2))*((2*P)-(2*(c/2)*(1-cos(th)))).*cos(th);            % Integral of {(DZ/Dx).cos(theta)) w.r.t theta
Int_4 = integral(fun4,theta,pi);                                              % limit is from theta to pi
Int_5 = Int_3 + Int_4;                                                        % Integral with limit 0 to pi
A1 = (2/pi)*(Int_5)                                                           % A1 formula

fun6 = @(th) (M/(P^2))*((2*P)-(2*(c/2)*(1-cos(th)))).*cos(2*th);              %Integral of {(DZ/Dx).cos(2*theta)) w.r.t theta
Int_6 = integral(fun6,0,theta);                                               %limit is 0 to theta

fun7 = @(th) (M/((1-P)^2))*((2*P)-(2*(c/2)*(1-cos(th)))).*cos(2*th);          %Integral of {(DZ/Dx).cos(2*theta)) w.r.t theta
Int_7 = integral(fun7,theta,pi);                                              % limit is from theta to pi
Int_8 = Int_6 + Int_7;                                                        % Integral limit from 0 to pi

A2 = (2/pi)*(Int_8)
CL = 2*pi*(A0 + (A1/2))                                                       % coefficient of lift
CM = (-pi/2)*(A0 + A1 - (A2/2))                                               % Coeffient of Moment
CP = (c/4)*(((2*A0) + (2*A1) -A2)/((2*A0)+ A1))                               % Center of pressure
AOA_L0 = ((Int/pi)-(Int_5/pi))                                                % Angle of attack alpha where Lift= 0 
AOA_L0_deg = AOA_L0*(180/pi)                                                  % in degrees


%% PLOTs


% Plotting figure of source panels

    figure(1);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    plot(XB,YB,'k-','LineWidth',3);                                         % Plot airfoil panels
    pB = plot(XB,YB,'ko','MarkerFaceColor','k');                            % Plot boundary points (black circles)
    pC = plot(XC,YC,'ko','MarkerFaceColor','m');                            % Plot control points (red circles)
    legend([pB,pC],{'Boundary Points','Control Points'});
	  xlabel('X ');                                                           % Set X-label
    ylabel('Y ');                                                           % Set Y-label
    xlim('auto');                                                           % Set X-axis limits to auto
    ylim('auto');                                                           % Set Y-axis limits to auto
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom


% Cp vectors at airfoil control points


    figure(2);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    Cps = abs(Cp*0.25);                                                     % Scale and make positive all Cp values
    for i = 1:1:length(Cps)                                                 % Loop over all panels
        X(1) = XC(i);                                                       % Control point X-coordinate
        X(2) = XC(i) + Cps(i)*cosd(betaD(i)+AoA);                           % Ending X-value based on Cp magnitude
        Y(1) = YC(i);                                                       % Control point Y-coordinate
        Y(2) = YC(i) + Cps(i)*sind(betaD(i)+AoA);                           % Ending Y-value based on Cp magnitude
        
        if (Cp(i) < 0)                                                      % If pressure coefficient is negative
            p{1} = plot(X,Y,'g-','LineWidth',2);                            % Plot 
        elseif (Cp(i) >= 0)                                                 % If pressure coefficient is zero or positive
            p{2} = plot(X,Y,'b-','LineWidth',2);                            % Plot 
        end
    end
    fill(XB,YB,'k');                                                        % Plot the airfoil as black polygon
	  legend([p{1},p{2}],{'Negetive Cp','Positive Cp'});                          % Show legend
    xlabel('X distance');                                                   % Set X-label
    ylabel('Cp');                                                           % Set Y-label
    xlim([-0.5 1.5]);                                                       % Set X-axis limits to auto
    ylim([-0.5 0.5]);                                                       % Set Y-axis limits to auto
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom

% Cp at upper and lower surface


    figure(3);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    midIndX = floor(length(xFoilCP)/2);                                     % Airfoil middle index for XFOIL data
    midIndS = floor(length(Cp)/2);                                          % Airfoil middle index for SPM data

    pSu = plot(XC(1:midIndS),Cp(1:midIndS),'ks','MarkerFaceColor','m');     % Plot Cp for upper surface of airfoil from SPM
    pSl = plot(XC(midIndS+1:end),Cp(midIndS+1:end),'ks',...                 % Plot Cp for lower surface of airfoil from SPM
                    'MarkerFaceColor','b');
    legend([pSu,pSl], {'Upper Surface','Lower Surface'});
    xlabel('X');                                                            % Set X-label
    ylabel('Cp');                                                           % Set Y-label
    xlim([0 1]);                                                            % Set X-axis limits
    ylim([-3 2]);                                                           % Set Y-axis limits to auto
    set(gca,'Ydir','reverse')                                               % Reverse direction of Y-axis
    title(['Airfoil: NACA 6412']);                                          % Set title
                                                          

