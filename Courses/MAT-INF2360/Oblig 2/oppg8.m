clear all;
close all;
clc;

for m = 1:4
    % spiller av bidrag fra resolusjonsrommet V0 for for b�de Haar-waveleten 
    % og v�r alternative stykkevise line�re wavelet som bruker phi som 
    % skaleringsfunksjon og psihat som moderwavelet
    %playDWTall(m);
    
    % spiller av detaljdelen for b�de Haar-waveleten og v�r alternative stykkevise
    % line�re wavelet
    playDWTalldifference(m);
end