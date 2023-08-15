%% Calculations
dt = 400;
G = 6.6738e-11;
AU = 149597870700;
m1 = 2e30;
m2 = 2*m1;
r1 = [2/3*AU 0];
r2 = [-1/3*AU 0];
v1 = [0 20000];
v2 = -1/2*v1;
i = 1;
number_of_calculations = 5*10^4;
while i <= number_of_calculations
    r = r2(i,:) - r1(i,:);
    F1 = G*m1*m2*(r)/norm(r)^3;
    F2 = -F1;
    v1 = v1 + (F1/m1)*dt;
    v2 = v2 + (F2/m2)*dt;
    r1(i+1,:) = r1(i,:) + v1*dt;
    r2(i+1,:) = r2(i,:) + v2*dt;
    i = i + 1;
end
%% Plots
s = size(r1);
di = 100;
i = 1;
r12 = m1*(m1+m2)/(m1*m2)*r1;
r21 = -r12;
while i <= s(1)
    hold off
    plot(r1(:,1),r1(:,2),'b:', r2(:,1),r2(:,2),'g:')
    axis([-2e11 2e11 -2e11 2e11])
    title('')
    hold on
    plot(r1(i,1),r1(i,2),'bo', r2(i,1),r2(i,2),'go')
    plot(0,0,'ro') % center of mass
    %plot([r1(i,1), r2(i,1)], [r1(i,2), r2(i,2)], 'r')
    
    r12_moved = r12(:,1) + r2(i,1); % moving the ellipses
    r12_moved(:,2) = r12(:,2) + r2(i,2);
    r21_moved = r21(:,1) + r1(i,1);
    r21_moved(:,2) = r21(:,2) + r1(i,2);
    plot(r12_moved(:,1), r12_moved(:,2), 'b:')
    plot(r21_moved(:,1), r21_moved(:,2), 'g:')
    
    i = i + di;
    pause(0.01)
end