function playDWTalldifference(m)

disp('Haar wavelet');
playDWTlowerdifference(m);

disp('Wavelet for piecewise linear functions, alternative version');
playDWTfilterslowerdifference(m,...
(1./128)*[-5 20 -1 -96 70 280 70 -96 -1 20 -5],...
(1./16)*[1 -4 6 -4 1],...
(1./16)*[1 4 6 4 1],...
(1./128)*[5 20 1 -96 -70 280 -70 -96 1 20 5]);

end