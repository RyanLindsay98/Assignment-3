clear
clc
%Ryan Lindsay 101038101
%% Constants %%

T =300;                  %Temp in K
K =1.38e-23;             %Boltsmann constant
Tmn =0.2e-12;            %mean time between collisions
Mo =9.11e-31;            %rest mass
Mn =0.26*Mo;             %effective mass of electrons
L =200e-09;               %Length of region
W =100e-09;               %Width of region
Pop =30000;                %number of particles
Vth = sqrt((K*T)/(Mn));   %Thermal velocity    
Tstep = 15e-15;           %time step of 15ns
q = 1.602e-19;             % electron charge
lengthE = 100;


%% Electron Modelling 

MFP = Tmn *Vth;

Ang = rand(Pop,1)*2*pi;  % Defines a random angle 

Pos = [rand(Pop,1)*L rand(Pop,1)*W Vth*cos(Ang) Vth*sin(Ang)];  %Creates an Array of particles with random X & Y positions and velocities 

initX = Pos(:,1); %The Initial X positions 

initY = Pos(:,2); % The initial Y positions 

colour = rand(Pop,1); % Ensures each electron will have its own colour 


% NEW TO Assignment 3:

Vy = 0;
Vx = 0.1;

EfieldX = Vx/L;
EfieldY = Vy/W;

EforceX = EfieldX*q;
EforceY = EfieldY*q;


Ax = zeros(Pop,1);
Ay = zeros(Pop,1);

Ay(:,1) = EforceY/Mn;
Ax(:,1) = EforceX/Mn;


cVec = zeros(lengthE,1);
eConc = 10^15;




initT = T; 

Pscat = 1- exp(-Tstep/Tmn);



Velo = sqrt(sum(Pos(:,3).^2)/Pop + sum(Pos(:,4).^2)/Pop);

probV = makedist('Normal', 'mu', 0, 'sigma', sqrt(K*T/Mn));

sumT = 0;



for i = 1 : lengthE      % Main Loop of the Function 
    
    % Probability of scattering 
    
    P = rand(Pop,1) < Pscat;
    
    Pos(P,3:4) = random(probV, [sum(P),2]);
    
    
    
    % Adding Acceleration in X and Y directions
    if Ax ~= 0
        
    
    Pos(:,3) = Pos(:,3) + Ax*Tstep;
    
    end
    
    
    
    if Ay ~=0
        
       
        Pos(:,4) = Pos(:,4) + Ay*Tstep;
        
    end
    
    
    
    
    Vrms = sqrt((Pos(:,3) .^ 2) + (Pos(:,4) .^ 2));
    
    Emob = mean(Vrms);
    
    
    driftV = Emob*EfieldX;
    
    avgCur = eConc*L*W*((driftV)/Pop)*q;
    
    cVec(i) = avgCur;
    
    
    
    
    % Finding velocity and Temperature 
    Velo = sqrt(sum(Pos(:,3).^2)/Pop + sum(Pos(:,4).^2)/Pop);
    
    newT = T + ((Mn * (Velo.^2) )/K/Pop/2);
    
    
    %Avg Temp Calc
    sumT = sumT + newT;
    
    avgT = sumT/i;
    
    
    
    newX = initX + Pos(:,3)*Tstep;    % The next X position of the particle
    
    newY = initY + Pos(:,4)*Tstep;     % The next Y position of the particle
    
    
    
    
    
    % Checking for Top and Bottom bounds 
    
    Yhigh = newY > W;
    newY(Yhigh) = 2*W - newY(Yhigh); 
    Pos(Yhigh,4) = -Pos(Yhigh,4);
    
    Ylow = newY < 0;
    newY(Ylow) = -newY(Ylow);
    Pos(Ylow,4) = -Pos(Ylow,4);
    
    
    
    % Checking for Left and Right Bounds 
    
    
    Xright = newX > L;
    newX(Xright) = newX(Xright) -L;
    
    
    Xleft = newX < 0;
    newX(Xleft) = newX(Xleft) + L;
    
    
    
    
    
    avgV=mean(sqrt(Pos(:,3).^2 + Pos(:,4).^2));
    MFP=(avgV*i)/1000000;
    
    
    % Plotting the movement of electrons
    
    
    Z = 500; % Step size of how many points are skipped
    
    sz = 10; % Size of bubbles
    
    
    figure(1)   
    
    scatter(initX(1:Z:end),initY(1:Z:end),sz,colour(1:Z:end),'filled');      
    hold on
    scatter(newX(1:Z:end), newY(1:Z:end),sz,colour(1:Z:end),'filled');         
    title('2D Trajectories')
    xlabel 'Length of Substrate'
    ylabel 'Width of Substrate' 
    axis([0 L 0 W]);
    
    
    
    
    
    
    %Re-initializing after 1 loop 
    initX = newX;
    initY = newY;
    initT = newT;
    
    
    

end

tVec = linspace(0,Tstep*lengthE,lengthE); 
%Plot of Current vs Time
    figure(3);
    plot( tVec,cVec)
    title('Average Current in X Direction ')
    xlabel('Time (s)')
    ylabel('Current (A)')
    hold on
 




% Density Plot 

n = 15; % Number of grid points 

[X Y] = meshgrid(0:L/n:L,0:W/n:W);  % Create a Grid 

Den = zeros(n+1,n+1);
Temp = zeros(n+1,n+1);

countT = 0;
countD = 0;

for i = 1:n
    
    
    XP1 = X(1,i);       
    XP2 = X(1,i+1);
    
    for j = 1:n
        
        YP1 = Y(j,1);
        YP2 = Y(j+1, 1);
        
        
          %check each frame
            for k=1:Pop
                
                if((Pos(k,1)>XP1 & Pos(k,1)<XP2) & Pos(k,2)<YP2 & Pos(k,2)>YP1)
                    
                    countT = countT+1;    
                    
                    Den(i,j) = Den(i,j)+1;      
                    
                    countD = countD + sqrt(Pos(k,3)^2+Pos(k,4)^2);  
                    
                    
                    if(countT >0)
                        Temp(i,j)=Mn*(countD^2)/(countT)/K/2;         
                    end
                end
            end
            
          countD=0;
           countT=0;
    end
        
    
end

%density map
    figure(5)
    surf(X,Y,Den)
    title 'Electron Density Map';
    zlabel 'Number of Electrons per Grid Point';
    ylabel 'Y Coordinate';
    xlabel 'X coordinate';
    
    
    
    %temperature map
    figure(6)
    surf(X,Y,Temp)
    title 'Temperature Density Map';
    zlabel 'Temperature per Grid Point';
    ylabel 'Y Coordinate';
    xlabel 'X coordinate';















