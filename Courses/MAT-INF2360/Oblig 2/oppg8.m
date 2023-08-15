clear all;
close all;
clc;

for m = 1:4
    % spiller av bidrag fra resolusjonsrommet V0 for for både Haar-waveleten 
    % og vår alternative stykkevise lineære wavelet som bruker phi som 
    % skaleringsfunksjon og psihat som moderwavelet
    %playDWTall(m);
    
    % spiller av detaljdelen for både Haar-waveleten og vår alternative stykkevise
    % lineære wavelet
    playDWTalldifference(m);
end