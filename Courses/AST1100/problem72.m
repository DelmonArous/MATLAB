% Constants
m1 = 6.4*10^23;         % mass of Mars [kg]
m2 = 100.0;             % mass of Beagle 2 [kg]           
radius = 3400.0*10^3;   % radius of Mars [m]
G = 6.67*10^-11;        % gravitational constant [N*(m/kg)^2]
k = 0.00016;            % friction constant [kg/s]            
n = 10^5;               % number of calculations
dt = 1;                 % timestep [s]

% Declaring arrays with initial values
r1 = zeros(n,2);
r2 = zeros(n,2);
v1 = zeros(n,2);
v2 = zeros(n,2);

r1(1,:) = [0 0];
r2(1,:) = [(-298-3400)*1000 0];
v1(1,:) = [0 0];
v2(1,:) = [0 -4000];

% Calculations using Euler-Cromer
i = 1; % counter
while i <= n
    rr = norm(r2(i,:)-r1(i,:)); % distance between the bodies 
    
    % acc. on Mars using N2L: necessary to find the next value 
    % of Mars' posisjon using Euler-Cromer so we can calculate 
    % the distance between the bodies once more
    a1 = (G*m2/rr^3)*r1(i,:);
    a2 = -(G*m1/rr^3)*r2(i,:)-(k/m2)*v2(i,:); % acc. on Beagle 2 using N2L
    
    v1(i+1,:) = v1(i,:) + dt*a1;
    v2(i+1,:) = v2(i,:) + dt*a2;
    
    r1(i+1,:) = r1(i,:) + dt*v1(i+1,:);
    r2(i+1,:) = r2(i,:) + dt*v2(i+1,:);
    
    if rr <= radius % Beagle 2 has landed
        break
    end
    i = i + 1;
    
end

% Surface approximation of Mars
theta = linspace(0,2*pi,n);
x = radius*cos(theta);
y = radius*sin(theta);
% Plot commands
hold on
plot(x,y)
plot(r2(:,1),r2(:,2),'--')
legend('Mars','Beagle 2 orbit')
xlabel('x-axis [m]')
ylabel('y-axis [m]')
axis('equal')
title('Trajectory of Beagle 2')
hold off