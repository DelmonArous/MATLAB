load -ascii motion1.d
load -ascii motion2.d

t1 = motion1(:,1);
n1 = length(t1);
dt1 = t1(1) - t1(2);
a1 = zeros(n1,2);
v1 = zeros(n1,2);
r1 = zeros(n1,2);
a1(:,1) = motion1(:,2);
a1(:,2) = motion1(:,3);
r1(1,:) = [0.0 0.0];
v1(1,:) = [0.0 0.0];

t2 = motion2(:,1);
n2 = length(t2);
dt2 = t2(1) - t2(2);
a2 = zeros(n2,2);
v2 = zeros(n2,2);
r2 = zeros(n2,2);
a2(:,1) = motion2(:,2);
a2(:,2) = motion2(:,3);
r2(1,:) = [0.0 0.0];
v2(1,:) = [0.0 0.0];

for i = 1:n1-1
    v1(i+1,:) = v1(i,:) + dt1*a1(i,:);
    r1(i+1,:) = r1(i,:) + dt1*v1(i,:);
end

for i = 1:n2-1
    v2(i+1,:) = v2(i,:) + dt2*a2(i,:);
    r2(i+1,:) = r2(i,:) + dt2*v2(i,:);
end

count1 = 0;
count2 = 0;
max_a1 = 0;
max_a2 = 0;

for i = 1:n1
    aks_mag = sqrt(a1(i,1)^2 + a1(i,2)^2); 
    if aks_mag > max_a1
        max_a1 = aks_mag;
        count1 = i;
    end
end

for i = 1:n2
    aks_mag = sqrt(a2(i,1)^2 + a2(i,2)^2); 
    if aks_mag > max_a2
        max_a2 = aks_mag;
        count2 = i;
    end
end

hold on
plot(r1(:,1), r1(:,2), '--')
plot(r2(:,1), r2(:,2), '-')
legend('Data-set: motion1.d', 'Data-set: motion2.d')
text(r1(count1,1), r1(count1,2),'\leftarrow max accel.','HorizontalAlignment','left')
text(r2(count2,1), r2(count2,2),'\leftarrow max accel.','HorizontalAlignment','left')
hold off
xlabel('x [m]');
ylabel('y [m]');
title('Plot of the data-set motion1.d and motion2.d')

% The motion from the data-set motion1.d indicates motion of an object 
% passing through a loop. In physical terms this may be a roller coaster 
% following a track which is leading to a vertical loop. In this case, 
% we have studied the motion of the cars attached to eachother as one 
% system. This because of the acceleration of each individual car 
% are equal.

% The motion from the data-set motion2.d indicates motion of an object 
% passing through a loop, but fails due to low acceleration at the
% beginning of the loop. The object makes it almost halfway through the
% loop before falling down and bounces a few times till the motion ends.
% The maximum acceleration is a result of the first bounce after the fall
% of the object.