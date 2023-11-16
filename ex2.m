%%%%Exercise 2

%%A
clear all;
clc;
close all;
pkg load signal;
pkg load statistics;
%%A1

%%A1
T = 10^-2;
over = 10;
Ts = T/over;
A = 4;
a = 0.5;

% make the phi function
[phi,t] = srrc_pulse(T, over, A, a);

Nf = 2048; % number of samples for the fourier signal

% making the F axis
Fs = 1/Ts;

% Nf samples of fourier axis 
F_axis = [-Fs/2 : Fs/Nf : Fs/2 - Fs/Nf ];

% Fourier transform of phi centralized at 0
XF = fftshift(fft(phi,Nf))*Ts;

% Plot the energy spectral density of XF using logarithmic y axis
figure(1)
semilogy(F_axis, abs(XF).^2);
xlim([-500, 500])
xlabel('F(Hz)')
title('Energy spectral density of srrc pulse');

%A2

N=100;

b = (sign (randn(N,1)) + 1) / 2 ; % bits eisodou.

X = bitsTo2PAM(b); % ta sumbola eisodou.


X_delta=(1/Ts) * upsample(X, over);

t_X=0:Ts:N*T-Ts;

figure(14);
stem(t_X,X_delta);
title('random symbols X');
%X(t)

Xt=conv(X_delta,phi)*Ts;
t_conv=min(t)+min(t_X) : Ts : max(t)+max(t_X);

figure(15);
plot(t_conv,Xt);
title('Convolution of Xd and phi signal');

%A3
Ttotal = t_conv(end) - t_conv(1); 

XFt = fftshift(fft(Xt,Nf)) * Ts;
numerator = abs(XFt.^2);
PxF = numerator / Ttotal;

figure(2);

plot(F_axis , PxF);
title('Periodgram of X(t) using plot');

figure(3);
semilogy(F_axis , PxF);
title('Periodgram of X(t) using semilogy');

K=500;
for i = 1 : K 
  b = (sign (randn(N,1)) + 1) / 2 ; % bits eisodou.

  X = bitsTo2PAM(b); % ta sumbola eisodou.

  X_delta=(1/Ts) * upsample(X, over);

  t_X=0:Ts:N*T-Ts;

  %X(t)

  Xt=conv(X_delta,phi)*Ts;
  t_conv=min(t)+min(t_X) : Ts : max(t)+max(t_X);

  Ttotal = t_conv(end) - t_conv(1); 
  
  XFt = fftshift(fft(Xt,Nf)) * Ts;
  numerator = abs(XFt.^2);
  PxFs(i,:) =   numerator / Ttotal; 
  
end

%%Plot some of the implementations of Xt , so we can see how the SxF looks like
figure(4);
subplot(2,2,1);
semilogy(F_axis , PxFs(200,:));
title('Periodgram of X(t) (for one implementation)');

subplot(2,2,2);
semilogy(F_axis , PxFs(300,:));
title('Periodgram of X(t) (for one implementation)');

subplot(2,2,3);
semilogy(F_axis , PxFs(400,:));
title('Periodgram of X(t) (for one implementation)');

subplot(2,2,4);
semilogy(F_axis , PxFs(500,:));
title('Periodgram of X(t) (for one implementation)');


%%Find the experimental SxF(psd) (arithmetic mean)

SxFexp = 1/500 * sum(PxFs); %%Sum all the columns of PxFs

%%Find the theoretical SxF = sigma^2/(t) * |F(Phi(f))|^2

SxFth = ( var(X)^2 ) / T * abs(XF).^2 ;

figure(5);

semilogy(F_axis , SxFexp ,F_axis ,  SxFth); 
title('2PAM');
legend('SxF experimental (N=100 , K=500 )','SxF theoretical');

%%%Do the same procedure , while increasing K,N


K=1000;
N= 200;

for i = 1 : K 
  b = (sign (randn(N,1)) + 1) / 2 ; % bits eisodou.

  X = bitsTo2PAM(b); % ta sumbola eisodou.

  X_delta=(1/Ts) * upsample(X, over);

  t_X=0:Ts:N*T-Ts;

  %X(t)

  Xt=conv(X_delta,phi)*Ts;
  t_conv=min(t)+min(t_X) : Ts : max(t)+max(t_X);

  Ttotal = t_conv(end) - t_conv(1); 
  
  XFt = fftshift(fft(Xt,Nf)) * Ts;
  numerator = abs(XFt.^2);
  PxFs2(i,:) =   numerator / Ttotal; 
  
end

%%Find the experimental SxF(psd) (arithmetic mean)

SxFexp2 = 1/1000 * sum(PxFs2); %%Sum all the columns of PxFs

figure(6);
semilogy(F_axis , SxFexp2 ,F_axis ,  SxFth); 
title('2PAM');
legend('SxF experimental(N=200 , K=1000 )','SxF theoretical');


%%%4 PAM
K=1000;
N= 200;

b1 = (sign (randn(N,1)) + 1) / 2 ; % bits eisodou.

X1 = bitsTo4PAM(b); % ta sumbola eisodou.

X_delta1=(1/Ts) * upsample(X1, over);

t_X=0:Ts:(N/2)*T-Ts;

%X(t)

Xt1=conv(X_delta1,phi)*Ts;
t_conv1=min(t)+min(t_X) : Ts : max(t)+max(t_X);
figure(16);
plot(t_conv1,Xt1);
title('Convolution of Xd and phi signal');



for i = 1 : K 
  b = (sign (randn(N,1)) + 1) / 2 ; % bits eisodou.

  X = bitsTo4PAM(b); % ta sumbola eisodou.

  X_delta=(1/Ts) * upsample(X, over);

  t_X=0:Ts:N*T-Ts;

  %X(t)

  Xt=conv(X_delta,phi)*Ts;
  t_conv=min(t)+min(t_X) : Ts : max(t)+max(t_X);

  Ttotal = t_conv(end) - t_conv(1); 
  
  XFt = fftshift(fft(Xt,Nf)) * Ts;
  numerator = abs(XFt.^2);
  PxFs4pam(i,:) =   numerator / Ttotal; 
  
end

SxFexp4pam = 1/1000 * sum(PxFs4pam); %%Sum all the columns of PxFs

%%Theoretical SxF 
SxFth1 = ( var(X1)^2 ) / T * abs(XF).^2 ;

figure(7);
semilogy(F_axis , SxFexp4pam ,F_axis ,  SxFth1); 
title('4PAM');
legend('SxF4pam experimental(N=200 , K=1000 )','SxF theoretical');


%%A5

%%The new phi signal will have 2*T and 2*over

[fi,t]=srrc_pulse(2*T, over*2, A, a);
%FI(F)
Ts = 2*T / (2*over) ; %% == Ts
FI=fftshift(fft(fi,Nf))*Ts;
N=100;
b = (sign(randn(N, 1)) + 1)/2;
Xn=bitsTo2PAM(b);
X_delta=(1/Ts) * upsample(Xn, 2*over);
t_X=0:Ts: N*2*T-Ts;

%X(t)

Xt=conv(X_delta,fi)*Ts;
t_conv=min(t)+min(t_X):Ts:max(t)+max(t_X);

Ttotal = t_conv(end) - t_conv(1); 

XFt = fftshift(fft(Xt,Nf)) * Ts;
numerator = abs(XFt.^2);
PxF = numerator / Ttotal;

figure(8);

plot(F_axis , PxF);
title('Periodgram of X(t) (Tnew = 2T) using plot');

figure(9);
semilogy(F_axis , PxF);
title('Periodgram of X(t) (Tnew = 2T) using semilogy');

K=500;
N=100;
for i = 1 : K 
  b = (sign (randn(N,1)) + 1) / 2 ; % bits eisodou.

  X = bitsTo2PAM(b); % ta sumbola eisodou.

  X_delta=(1/Ts) * upsample(X, 2*over);

  t_X=0:Ts:N*T-Ts;

  %X(t)

  Xt=conv(X_delta,phi)*Ts;
  t_conv=min(t)+min(t_X) : Ts : max(t)+max(t_X);

  Ttotal = t_conv(end) - t_conv(1); 
  
  XFt = fftshift(fft(Xt,Nf)) * Ts;
  numerator = abs(XFt.^2);
  PxFs3(i,:) =   numerator / Ttotal; 
  
end

%%Plot some of the implementations of Xt , so we can see how the SxF looks like
figure(10);
subplot(2,2,1);
semilogy(F_axis , PxFs3(200,:));
title('Periodgram of X(t) (Tnew = 2T) (for one implementation)');

subplot(2,2,2);
semilogy(F_axis , PxFs3(300,:));
title('Periodgram of X(t) (Tnew = 2T) (for one implementation)');

subplot(2,2,3);
semilogy(F_axis , PxFs3(400,:));
title('Periodgram of X(t) (Tnew = 2T) (for one implementation)');

subplot(2,2,4);
semilogy(F_axis , PxFs3(500,:));
title('Periodgram of X(t) (Tnew = 2T) (for one implementation)');


%%Find the experimental SxF(psd) (arithmetic mean)

SxFexp3 = 1/500 * sum(PxFs3); %%Sum all the columns of PxFs

%%Find the theoretical SxF = sigma^2/(t) * |F(Phi(f))|^2

SxFth = ( var(X)^2 ) / (2*T) * abs(XF).^2 ;

figure(11);

semilogy(F_axis , SxFexp3 ,F_axis ,  SxFth); 
title('2PAM (Tnew = 2T)');
legend('SxF experimental (N=100 , K=500 )','SxF theoretical');

%%%Do the same procedure , while increasing K,N


K=1000;
N= 200;

for i = 1 : K 
  b = (sign (randn(N,1)) + 1) / 2 ; % bits eisodou.

  X = bitsTo2PAM(b); % ta sumbola eisodou.

  X_delta=(1/Ts) * upsample(X, 2*over);

  t_X=0:Ts:N*T-Ts;

  %X(t)

  Xt=conv(X_delta,phi)*Ts;
  t_conv=min(t)+min(t_X) : Ts : max(t)+max(t_X);

  Ttotal = t_conv(end) - t_conv(1); 
  
  XFt = fftshift(fft(Xt,Nf)) * Ts;
  numerator = abs(XFt.^2);
  PxFs4(i,:) =   numerator / Ttotal; 
  
end

%%Find the experimental SxF(psd) (arithmetic mean)

SxFexp4 = 1/1000 * sum(PxFs4); %%Sum all the columns of PxFs

figure(12);
semilogy(F_axis , SxFexp4 ,F_axis ,  SxFth); 
title('2PAM (Tnew = 2T)');
legend('SxF experimental(N=200 , K=1000 )','SxF theoretical');


%%%%%B

%B1

%%Y(t) = X cos(2pifot + theta), X~N(0,1) , theta ~ U [0,2pi)
f0 = 100 ;
t = [0: 0.00001 : 0.05];

%%plot 5 implementations 

figure();

for i=1 : 5 
  x = normrnd( 0 , 1 ); 
  xs(i) = x;
  
  theta = unifrnd(0,2*pi);
  thetas(i) = theta;
  
  y = x * cos(2*pi*f0*t + theta);
  plot(t,y);
  title('5 implementations of Y(t) = X cos(2*pi*fo*t + theta)');
  hold on;
end 
str1 = sprintf('Implementation of Y(t) with x= %.2f and theta= %.2f', xs(1) , thetas(1)) ;
str2 = sprintf('Implementation of Y(t) with x= %.2f and theta= %.2f', xs(2) , thetas(2)) ;
str3 = sprintf('Implementation of Y(t) with x= %.2f and theta= %.2f', xs(3) , thetas(3)) ;
str4 = sprintf('Implementation of Y(t) with x= %.2f and theta= %.2f', xs(4) , thetas(4)) ;
str5 = sprintf('Implementation of Y(t) with x= %.2f and theta= %.2f', xs(5) , thetas(5)) ;
legend(str1,str2,str3,str4,str5);
hold off;





 