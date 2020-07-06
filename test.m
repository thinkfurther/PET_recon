clf;clear all;
delta_t = 0.001;
t = 0 : delta_t : 5 - delta_t;
N = length(t);

w0 = 2*pi*1;
y = cos(w0*t);
delta_w = 1 / delta_t;
w = delta_w * (-N/2 : N/2 -1) /N;

y_fft = fftshift(fft(y))/N;
figure(1);
subplot(1,2,1);
plot(t, y);
subplot(1,2,2);
plot(w, abs(y_fft));