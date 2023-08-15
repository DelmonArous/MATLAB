antsec = 4/400;
fs = 44100;
t1 = linspace(0, 4./440, fs*(4./440));
t2 = linspace(4/440, 12./440, fs*(8./440));
t3 = linspace(12/440, 20./440, fs*(8./440));

f1 = 0*t1;
f2 = 2*((440*t2-4)/8).*sin(2*pi*440*t2);
f3 = 2*sin(2*pi*440*t3);

S = [f1 f2 f3];
S = S/max(abs(S));
playerobj = audioplayer(S,fs);
playblocking(playerobj);

% samplesperperiod = round(200*fs/440);
% oneperiod = zeros(1,samplesperperiod);
% t = linspace(0, 200/440, samplesperperiod);
% 
% for i = 1:samplesperperiod
%     if (0 <= i && i < 40);
%         oneperiod(i) = 0;
%     elseif (40 <= i && i < 120);
%         oneperiod(i) = (110*t(i) - 1).*sin(2*pi*440*t(i));
%     elseif (120 <= i && i <= 200);
%         oneperiod(i) = 2.*sin(2*pi*440*t(i));
%     end
% end
% 
% allsamples = zeros(1,antsec*2.2*length(oneperiod));
% 
% for k = 1:(antsec*2.2)
%     allsamples(((k-1)*length(oneperiod)+1):k*length(oneperiod))=oneperiod;
% end
% 
% playerobj = audioplayer(allsamples,fs);
% playblocking(playerobj);