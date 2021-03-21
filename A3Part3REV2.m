clear
clc
%% Constants %%

T =300;                  %Temp in K
K =1.38e-23;             %Boltsmann constant
Tmn =0.2e-12;            %mean time between collisions
Mo =9.11e-31;            %rest mass
Mn =0.26*Mo;             %effective mass of electrons
L =200e-09;               %Length of region
W =100e-09;               %Width of region
Pop =100000;                %number of particles
Vth = sqrt((K*T)/(Mn));   %Thermal velocity    
Tstep = 15e-15;           %time step of 15ns
q = 1.602e-19;             % electron charge
lengthE = 50;


%% Electron Modelling 

MFP = Tmn *Vth;

Ang = rand(Pop,1)*2*pi;  % Defines a random angle 

Pos = [rand(Pop,1)*L rand(Pop,1)*W Vth*cos(Ang) Vth*sin(Ang)];  %Creates an Array of particles with random X & Y positions and velocities 

initX = Pos(:,1); %The Initial X positions 

initY = Pos(:,2); % The initial Y positions 

colour = rand(Pop,1); % Ensures each electron will have its own colour 



%Box Definition 

xBox = [L/3 2*L/3];

yBox1 = [5*W/8 W];

yBox2 = [0 3*W/8];





% Checking For Particles inside Boxes

%Top Box

 bCol = (initX(:,1)> xBox(1) & initX(:,1) < xBox(2) & initY(:,1) > yBox1(1) & initY(:,1) < yBox1(2));
 
 
while(sum(bCol) > 0)
    
initX(bCol) = rand(sum(bCol),1)*L;
initY(bCol) = rand(sum(bCol),1)*W;


bCol = (initX(:,1)> xBox(1) & initX(:,1) < xBox(2) & initY(:,1) > yBox1(1) & initY(:,1) < yBox1(2));


end


% Bottoma Box

bCol = (initX(:,1)> xBox(1) & initX(:,1) < xBox(2) & initY(:,1) > yBox2(1) & initY(:,1) < yBox2(2));

while(sum(bCol) > 0)
    
initX(bCol) = rand(sum(bCol),1)*L;
initY(bCol) = rand(sum(bCol),1)*W;


bCol = (initX(:,1)> xBox(1) & initX(:,1) < xBox(2) & initY(:,1) > yBox2(1) & initY(:,1) < yBox2(2));


end





% NEW TO Assignment 3:

Vy = 0;
Vx = 0.8;
EfieldX = Vx/L;
EfieldY = Vy/W;
Ay = zeros(Pop,1);
Ax = zeros(Pop,1);
Ay = zeros(Pop,1);
Ay = EfieldY*q/Mn;
Ax(:,1) = EfieldX*q/Mn;


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
    
    xBox = [L/3 2*L/3];

    yBox1 = [5*W/8 W];

    yBox2 = [0 3*W/8];

    
    if i > 1
    
    % Checking For Box Values 
    
        bCol = (newX>xBox(1) & newX <xBox(2) & newY < yBox2(2)) | (newX>xBox(1) & newX <xBox(2) & newY > yBox1(1)); % Definition of the two boxes
    
    
    if sum(bCol) ~= 0       % if bCol contains a value, we have an electron at a barrier 
        
        if ((newX(bCol) > xBox(1)) & newX(bCol) < xBox(2));     % Checking X plane 
            Pos(bCol,1) = prevX(bCol,1);
            Pos(bCol,4) = -Pos(bCol,4);
            Pos(bCol,3)=Vth*cos(randn()*2*pi);
            
            
        else
            Pos(bCol,2) = prevY(bCol,2);
            Pos(bCol,3) = -Pos(bCol,3);    % Checking Y plane 
            
            Pos(bCol,4)=Vth*sin(randn()*2*pi);
            
        end
    end
    
    
    end
    
    
    
    
    % Plotting the movement of electrons
    
    
    Z = 1000; % Step size of how many points are skipped
    
    sz = 10; % Size of bubbles
    
%     
%     figure(1)   
%     
%     scatter(initX(1:Z:end),initY(1:Z:end),sz,colour(1:Z:end),'filled');      
%     hold on
%     scatter(newX(1:Z:end), newY(1:Z:end),sz,colour(1:Z:end),'filled');         
%     rectangle( 'Position', [xBox(1) 0 (xBox(2)-xBox(1)) yBox2(2)])
%     rectangle('Position', [xBox(1) yBox1(1) (xBox(2)-xBox(1)) yBox2(2)])
%     xlabel 'Length of Substrate'
%     ylabel 'Width of Substrate' 
%     axis([0 L 0 W]);
    
    
   
    
    
    
    %Re-initializing after 1 loop 
    initX = newX;
    initY = newY;
    initT = newT;
    
    
    prevX = newX;
    prevY = newY;
    
    
    xBox = [L/3 2*L/3];

    yBox1 = [5*W/8 W];

    yBox2 = [0 3*W/8];
    
    
    figure(7)
    Xelec = linspace(0, 200e-9, 1000);
    Yelec = linspace(0, 100e-9, 500);
    n = histcounts2(newY, newX, Yelec, Xelec);
    imagesc(Xelec,Yelec,n),colorbar,title('Electron Density Map');
    
    
    
    
    
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

n = 30; % Number of grid points 

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
%     figure(5)
%     surf(X,Y,Den)
%     title 'Electron Density Map';
%     zlabel 'Number of Electrons per Grid Point';
%     ylabel 'Y Coordinate';
%     xlabel 'X coordinate';
    
    
    
%    temperature map
%     figure(6)
%     surf(X,Y,Temp)
%     title 'Temperature Density Map';
%     zlabel 'Temperature per Grid Point';
%     ylabel 'Y Coordinate';
%     xlabel 'X coordinate';
    
    
    
    