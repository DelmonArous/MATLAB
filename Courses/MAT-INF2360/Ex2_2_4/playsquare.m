function playsquare(T)

antsec = 3;
fs = 44100;
samplesperperiod = round(fs*T);
oneperiod = [ones(1,round(samplesperperiod/2))...
    -ones(1,round(samplesperperiod/2))];
allsamples = zeros(1,(antsec/T)*length(oneperiod));

for k = 1:(antsec/T)
    allsamples(((k-1)*length(oneperiod)+1):k*length(oneperiod))=oneperiod;
end

playerobj = audioplayer(allsamples, fs);
playblocking(playerobj);

end