% Monte Carlo Modelling of Electron Transport

close all
clear 
format short
clc

% Constants
m0 = 9.11e-31; %kg
me = 0.26*m0;
q = 1.602e-19; %C
kB = 1.38066e-23; %J/K
nm = 1e-9; %nanometre
ps = 1e-12; %picosecond

% Dimensions
% Elec = 1; % simulates for 1 particle
Elec = 3; % simulates for 5 particles at once
xdim = 200; %nm
ydim = 100; %nm and need to make same length
x = zeros(1,xdim)*nm;
y = zeros(1,ydim*2)*nm; %need to make same length 
% Reg = zeros(xdim*1e9,ydim*1e9)*nm;  % semiconductor region

% Specifics
InitialTemp = 300; % in [K]
vtFirst = sqrt(2*kB*InitialTemp/me); %1.8701e5 in m/s
Tmn = 0.2*ps; 

% Simulation
steps = 10; 
% t = zeros(1,steps);
% Temp = zeros(1,steps);
dt = nm/vtFirst; %5.347e-15s
Temp = zeros(1,Elec);
t = zeros(1,Elec);

on = 1;
off = 0;
scatter = on; %scattering on/off switch
% scatter = off; %scattering on/off switch
collide = 0; % counter for collisons
p = 1 - exp(-dt/Tmn); % probability for scatter threshold
re = 0; % default rebound velocity factor
% re = -0.5; % Bounces back at half its initial speed

for i = 1:steps    
   
    if i == 1
        t(i,:) = dt;
        Temp(i,:) = InitialTemp;
        vt = vtFirst;        
        
        Boxes = {};
        
        Boxes{1}.X = [0.9 1.1]*100*nm;
        Boxes{1}.Y = [0.8 1.0]*100*nm;
%         Boxes{1}.BC = 0.0;
        
        Boxes{2}.X = [-.9 1.1]*100*nm;
        Boxes{2}.Y = [0.0 0.2 1.1]*100*nm;
%         Boxes{2}.BC = 0.0; s
        
        e=1;
%         while e <= Elec %random positions for each particle
                        
            randx = round(rand*xdim)*nm;
            randy = round(rand*ydim)*nm;
            
%             % check to see if in box
%             if Boxes{1}.X(1,1)<randx && randx <Boxes{1}.X(1,2) &&...
%                 Boxes{1}.Y(1,1)<randy && randy<Boxes{1}.Y(1,2) &&...
%                 Boxes{2}.X(1,1)<randx && randx<Boxes{2}.X(1,2) &&...
%                 Boxes{2}.Y(1,1)<randy && randy<Boxes{2}.Y(1,2) 
%                                 
%             else
                x(i,e) = randx;
                y(i,e) = randy;
                e=e+1;
%             end
                      
%         end 
        
    else
        t(i,:) = t(i-1,:) + dt;
        r = rand(1,Elec); % Chance of scattering
        
        % Particle direction
        if i == 2 % intial angle, regardless of scattering
            
            for e = 1:Elec
                theta = rand*360; % random angle in deg
                dz = vtFirst*Tmn; % straight path
                dx(i,e) = sind(theta)*dz; % new x dif
                dy(i,e) = cosd(theta)*dz; % new y dif 
                x(i,e) = x(i-1,e) + dx(i,e);
                y(i,e) = y(i-1,e) + dy(i,e); 
                Temp(i,e) = Temp(i-1,e);
                vtNew = vtFirst;
            end

        elseif scatter == off % no scattering, just continue along inital path
            
            for e = e:Elec
                dx(i,e) = x(2,e)-x(1,e);
                dy(i,e) = y(2,e)-y(1,e); 
            end
             
            vtNew = vt;
            Temp(i,:) = Temp(i-1,:);

        else % yes scatter
            
            for e = 1:Elec 
               
                if r(1,e) < p*100   % check each electron
                    
                    scale = sqrt(kB*Temp(i-1,e)/me); %scaling factor
                    dof = 3; %degrees of freedom
%                     vtNew = vtFirst*chi2rnd(dof);
                    vtNew(i,e) = vtFirst*(rand/100+0.995);

                    theta = rand*360; % random angle in deg
                    dz(i,e) = vtNew(i,e)*Tmn; % straight path
                    dx(i,e) = sind(theta)*dz(i,e); % new x difference
                    dy(i,e) = cosd(theta)*dz(i,e); % new y difference 
                    
                    % updating position, velocity, temperature
                    TempNew = me*vtNew.^2/(2*kB); 
                    Temp(i,e) = TempNew(i,e);                    
                    x(i,e) = x(i-1,e) + dx(i,e);
                    y(i,e) = y(i-1,e) + dy(i,e); 
                    vt = vtNew;
                    
                    collide = collide + 1;
                    
                end
                
            end
                                            
        end
       
                
        %Boundary Conditions
%         for e = 1:Elec
%             
%             if x(i,:) < 0
%             x(i,:) = x(i,:) + 2e-7;            
%             elseif x(i,:) > xdim*nm
%                 x(i,:) = x(i,:) - 2e-7;
%             end
% 
%             if y(i,:) < 0
%                 y(i,:) = x(i,:) + 2e-7;
%             elseif y(i,:) > ydim*nm
%                 y(i,:) = x(i,:) - 2e-7;
%             end
% 
%         end
                
    end

%     1c-i)  Particle Plot   
    for e = 1:Elec
        figure(1)
        plot(x(:,e),y(:,e));
        grid on;
        hold on; 
    end
    xlim([0 xdim*nm])
    ylim([0 ydim*nm])
    xlabel('X (m)')
    ylabel('Y (m)')
    title(['Time Passed t: ', num2str(t(i)/ps), ...
        'ps Collsions: ', num2str(collide)]) 

    % 1c-ii) Temp plot  
    for e = 1:Elec
        figure(2)
        plot(t(:,e),Temp(:,e));
        grid on;
        hold on;
    end 
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    title(['Current Temperature: ', num2str(mean(Temp(i,:))), ...
        'K Max T: ', num2str(max(Temp,[],'all')),'K Min T: ',...
        num2str(min(Temp,[],'all')),'K']) % change T
        
         
    pause(0.05)
end

display('Seconds Passed', num2str(t(i)));
finalTemp = num2str(mean(Temp(i,:)))



