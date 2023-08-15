b = 0.1;    % m
g = 9.8;    % m/s^2
time = 2;   % s
dt = 0.01;  % s
u = 0.1;    % m/s
my_s = 0.6; 
my_d = 0.3;

n = ceil(time/dt);

v0 = 0.1;   % m/s
x0 = 0.0;   % m

t = zeros(n,1);
x = zeros(n,1);
v = zeros(n,1);
Fx = zeros(n,1);

v(1) = v0;
x(1) = x0;

for m = [0.1 1.0]
    
    if m == 0.1;
        line = '-';
        k = 10;
    else
        line = '--';
        k = 100;
    end
    
    for i = 1:n-1
        x_b = u*t(i) + x0 + b;
        N = m*g;
        Fx(i) = k*(x_b - x(i) - b);
        if v(i) == 0
            f = -Fx(i);
            if abs(f) > my_s*N
                F = Fx(i) - sign(f)*my_d*N;
            else
                F = 0;
            end
        else
            F = Fx(i) - sign(v(i))*my_d*N;
        end
        a = F/m;
        v(i+1) = v(i) + a*dt;
        if (v(i) ~= 0.0) && (sign(v(i+1))~=sign(v(i)))
            v(i+1) = 0.0;
        end
        x(i+1) = x(i) + v(i+1)*dt;
        t(i+1) = t(i) + dt;
    end
    
    hold on
    plot(t, Fx, line)
    xlabel('t [s]')
    ylabel('F [N]')
    title('Plot of the spring force on the block')
    hold off
    
end

legend('m = 0.1kg, k = 10N/m','m = 1.0kg, k = 100N/m')
