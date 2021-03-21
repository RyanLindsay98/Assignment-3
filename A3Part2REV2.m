%% Ryan Lindsay 101038101
clear 
clc
x = 3;  %Dimensions 
y = 2;


dx = 0.1;
dy = 0.1;

nx = 200;
ny = 100;

% Create matrix 
G = sparse(nx*ny, nx*ny);
F = zeros(nx*ny,1);


%% Constants

T =300;                  %Temp in K
K =1.38e-23;             %Boltsmann constant
Tmn =0.2e-12;            %mean time between collisions
Mo =9.11e-31;            %rest mass
Mn =0.26*Mo;             %effective mass of electrons
L =200e-09;               %Length of region
W =100e-09;               %Width of region
Pop =1000;                %number of particles
Vth = sqrt((K*T)/(Mn));   %Thermal velocity    
Tstep = 5e-15;           %time step of 15ns
q = 1.602e-19;             % electron charge
lengthE = 1000;

%% Finite Differences 



cMap = zeros(nx,ny)+1;
V = zeros(nx,ny);
F = zeros(nx*ny,1);

for i = 1:nx
    for j =1:ny
        
        if j <= ny/3 || j>= 2*ny/3
            
            if i >= nx/3 && i<= 2*nx/3
                
                cMap(i,j) = 10^-2;
                
            end
        end
    end
end


for i = 1:nx
        for j = 1:ny
            n = j+(i-1)*ny;
            
            
            if (i == 1)
                F(n) = 1;
                G(n,n) = 1;
                
            elseif (i == nx)
                F(n) = 0;
                G(n,n) = 1;
                
            elseif (j == 1)
                
                npy = j+1 + (i-1)*ny;
                npx = j + (i-1+1)*ny;
                nmx = j + (i-1-1)*ny;
                nmy = j-1 + (i-1)*ny;

                
                rpx = (cMap(i,j) + cMap(i+1,j))/2;
                rmx = (cMap(i,j) + cMap(i-1,j))/2;
                rpy = (cMap(i,j) + cMap(i,j+1))/2;

                G(n,n) = -(rpx+rmx+rpy);
                G(n,npx) = rpx; 
                G(n,nmx) = rmx;
                G(n,npy) = rpy;
                
            elseif (j == ny)
                
                nmy = j-1 + (i-1)*ny;
                npx = j + (i-1+1)*ny;
                nmx = j + (i-1-1)*ny;

                
                rpx = (cMap(i,j) + cMap(i+1,j))/2;
                rmx = (cMap(i,j) + cMap(i-1,j))/2;
                rmy = (cMap(i,j) + cMap(i,j-1))/2;

                G(n,n) = -(rpx+rmx+rmy);
                G(n,npx) = rpx; 
                G(n,nmx) = rmx;
                G(n,nmy) = rmy;
            else
                
                npy = j+1 + (i-1)*ny;
                nmy = j-1 + (i-1)*ny;
                npx = j + (i-1+1)*ny;
                nmx = j + (i-1-1)*ny;

                
                rpx = (cMap(i,j) + cMap(i+1,j))/2;
                rmx = (cMap(i,j) + cMap(i-1,j))/2;
                rpy = (cMap(i,j) + cMap(i,j+1))/2;
                rmy = (cMap(i,j) + cMap(i,j-1))/2;

                
                G(n,n) = -(rpx+rmx+rpy+rmy);
                G(n,npx) = rpx; 
                G(n,nmx) = rmx;
                G(n,npy) = rpy;
                G(n,nmy) = rmy;
            end
        end
    end

    
    V = G\F;
    Vrs = reshape(V,ny,nx);
    figure(1)
    surf(Vrs)
    title('Voltage Plot')
    xlabel(' X Direction')
    ylabel(' Y Direction')
    zlabel('Votage (V)')

    
    [Ex Ey] = gradient(-Vrs);
    figure(2)
    quiver(Ex,Ey)
    title('Quiver Plot')
    xlabel(' X Direction')
    ylabel(' Y Direction')
    
    
    






%% Monte Carlo

MFP = Tmn *Vth;

Ang = rand(Pop,1)*2*pi;  % Defines a random angle 

Pos = [rand(Pop,1)*L rand(Pop,1)*W Vth*cos(Ang) Vth*sin(Ang)];  %Creates an Array of particles with random X & Y positions and velocities 

initX = Pos(:,1); %The Initial X positions 

initY = Pos(:,2); % The initial Y positions 

colour = rand(Pop,1); % Ensures each electron will have its own colour 





%Field in X direction
ForceX = Ex'*q;

Ax = 1e9*ForceX/Mn;



%Field in Y direction

ForceY = Ey'*q;

Ay = 1e9*ForceY/Mn;



% Concentration
cVec = zeros(lengthE,1);
eConc = 10^15;




initT = T; 

Pscat = 1- exp(-Tstep/Tmn);

Velo = sqrt(sum(Pos(:,3).^2)/Pop + sum(Pos(:,4).^2)/Pop);

probV = makedist('Normal', 'mu', 0, 'sigma', sqrt(K*T/Mn));

sumT = 0;



%Box Definition 

xBox = [L/3 2*L/3];

yBox1 = [2*W/3 W];

yBox2 = [0 W/3];





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



for i = 1 : lengthE      % Main Loop of the Function 
    
    % Probability of scattering 
    
    P = rand(Pop,1) < Pscat;
    
    Pos(P,3:4) = random(probV, [sum(P),2]);
    
    
    newX = initX + Pos(:,3)*Tstep;    % The next X position of the particle
    
    newY = initY + Pos(:,4)*Tstep;     % The next Y position of the particle
    
    
    dummyX(:,i) = Pos(:,3);
    
    dummyY(:,i) = Pos(:,4);
    
    
    [xB, Xedge] = discretize(dummyX(:,i),200);
    [yB, Yedge] = discretize(dummyY(:,i),100);
    
    
    dummyX(:,1) = dummyX(:,1) + (1/2)*Ax(sub2ind(size(Ax),xB,yB))*Tstep;
    dummyY(:,1) = dummyY(:,1) + (1/2)*Ay(sub2ind(size(Ay),xB,yB))*Tstep;
    
    
    
    Pos(:,3) = dummyX(:,i);
    Pos(:,4) = dummyY(:,i);
    
    
    newX = initX + Pos(:,3)*Tstep;
    newY = initY + Pos(:,4)*Tstep;
    
    

    
    
    
    
    % Finding velocity and Temperature 
    Velo = sqrt(sum(Pos(:,3).^2)/Pop + sum(Pos(:,4).^2)/Pop);
    
    newT = T + ((Mn * (Velo.^2) )/K/Pop/2);
    
    
    %Avg Temp Calc
    sumT = sumT + newT;
    
    avgT = sumT/i;
    
    
    
    
    
    
    
    
    
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
    
    
    
    
    
    avgV=mean(sqrt(Pos(:,3).^2 + Pos(:,4).^2));
    MFP=(avgV*i)/1000000;
    
    
    % Plotting the movement of electrons
    
    
    Z = 100; % Step size of how many points are skipped
    
    sz = 10; % Size of bubbles
    
    
    figure(6)   
    
    scatter(initX(1:Z:end),initY(1:Z:end),sz,colour(1:Z:end),'filled');      
    hold on
    scatter(newX(1:Z:end), newY(1:Z:end),sz,colour(1:Z:end),'filled');         
    rectangle( 'Position', [xBox(1) 0 (xBox(2)-xBox(1)) yBox2(2)])
    rectangle('Position', [xBox(1) yBox1(1) (xBox(2)-xBox(1)) yBox2(2)])
    xlabel 'Length of Substrate'
    ylabel 'Width of Substrate' 
    title('2D Trajectories')
    axis([0 L 0 W]);
    
    
    
    
    
    
    %Re-initializing after 1 loop 
    initX = newX;
    initY = newY;
    initT = newT;
    
    
     prevX = newX;
    prevY = newY;
    
    
    xBox = [L/3 2*L/3];

    yBox1 = [2*W/3 W];

    yBox2 = [0 W/3];
    
    
end


