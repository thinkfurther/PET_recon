clear all;
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

f = Fs*(-(L/2):(L/2 - 1))/L;

Y = fft(S);
Y = fftshift(Y);
P2 = abs(Y/L);

figure(2);
plot(f,P2) 
title('Double-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')