%%%%%
%
% Implementation of a Jansen-Rit model, both with the original sigmoidal 
% relationship between mean membrane potential and mean output spike rate
% and with a linearised relationship.
%
%%%%%

clear;

pkg load statistics

% Parameters, as defined in paper text
a = 100;
A = 0.00325;
b = 50;
B = 0.022;
v0 = 0.006;
r = 560;
e0 = 2.5;
C = 300;
C1 = C;
C2 = 0.8*C1;
C3 = C1/4;
C4 = C1/4;

% Input noise
PLwrLmt = 120;
PUpprLmt = 320;

% Duration of simulation
MxTm = 4;

% Time step
TmStp = 0.0001;


function [y0a, y0b, y1, y2, y3a, y3b, y4, y5] = JansenRit(y0a, y0b, y1, y2, y3a, y3b, y4, y5, C1, C2, C3, C4, a, A, b, B, e0, r, v0, P, TmStp, n, Lnrsd)

  dy0adt = y3a(n-1);
  dy0bdt = y3b(n-1);
  dy1dt = y4(n-1);
  dy2dt = y5(n-1);
  dy3adt = A*a*C1*sgmd(y1(n-1)-y2(n-1), v0, r, e0, Lnrsd)-2*a*y3a(n-1)-a^2*y0a(n-1);
  dy3bdt = A*a*C3*sgmd(y1(n-1)-y2(n-1), v0, r, e0, Lnrsd)-2*a*y3b(n-1)-a^2*y0b(n-1);
  dy4dt = A*a*(P(n)+C2*sgmd(y0a(n-1), v0, r, e0, Lnrsd))-2*a*y4(n-1)-a^2*y1(n-1);
  dy5dt = B*b*(C4*sgmd(y0b(n-1), v0, r, e0, Lnrsd))-2*b*y5(n-1)-b^2*y2(n-1);

  y0a(n) = y0a(n-1) + TmStp*dy0adt;
  y0b(n) = y0b(n-1) + TmStp*dy0bdt;
  y1(n) = y1(n-1) + TmStp*dy1dt;
  y2(n) = y2(n-1) + TmStp*dy2dt;
  y3a(n) = y3a(n-1) + TmStp*dy3adt;
  y3b(n) = y3b(n-1) + TmStp*dy3bdt;
  y4(n) = y4(n-1) + TmStp*dy4dt;
  y5(n) = y5(n-1) + TmStp*dy5dt;

end


function S = sgmd(v, v0, r, e0, Lnrsd)

  if Lnrsd == 0

    S = 2*e0/(1+exp(r*(v0-v))) - 2*e0/(1+exp(r*(v0)));

  elseif Lnrsd == 1

    S = v*(2*e0*r*exp(r*v0))/(1+exp(r*v0))^2;

  end

end


function [f, FT] = FrrTrnsfrm(Tm, yr)

    TmStp = Tm(2) - Tm(1);

    [vrbl0, NmbrPnts] = size(yr);

    FT = fft(yr);
    FT = fftshift(FT);
    f = (-(NmbrPnts-1)/2:(NmbrPnts-1)/2)*1/TmStp*1/NmbrPnts;
    FT = abs(FT/NmbrPnts);
    [vrbl0, IndxMn] = min(abs(f-1));
    [vrbl0, IndxMx] = min(abs(f-60));
    f = f(IndxMn:IndxMx);
    FT = FT(1,IndxMn:IndxMx);

end


function [OutputTF] = TrnsfrFnctn(Prmtrs, fInptFT)

  a = Prmtrs(1,1);
  A = Prmtrs(2,1);

  b = Prmtrs(3,1);
  B = Prmtrs(4,1);

  r = Prmtrs(5,1);
  e0 = Prmtrs(6,1);
  v0 = Prmtrs(7,1);

  C1 = Prmtrs(8,1);
  C2 = Prmtrs(9,1);
  C3 = Prmtrs(10,1);
  C4 = Prmtrs(11,1);

  gmm = (2*e0*r*exp(r*v0))/((1+exp(r*v0))^2);

  [vrbl0 NmbrPnts] = size(fInptFT);

  Frqncs = fInptFT;

  %%%
  %
  % Note the factor of 2pi introduced here to equate with FFT
  %
  %%%

  for n=1:NmbrPnts

      s = 1i * 2 * pi * Frqncs(n);

      Trm1 = ((s+a)^2)/(A*a);
      Trm2 = (a*A*C1*C2*gmm^2)/((s+a)^2);
      Trm3 = (b*B*C3*C4*gmm^2)/((s+b)^2);

      TF(n) = 1/(Trm1 - Trm2 + Trm3);
      OutputTF(n) = TF(n);
      OutputTF(n) = (abs(OutputTF(n)));
      fTF(n) = fInptFT(n);

  end

end


FntSz = 15;
figure(1);

% Sigmoid function

y0a(1:2) = 0;
y0b(1:2) = 0;
y1(1:2) = 0;
y2(1:2) = 0;
y3a(1:2) = 0;
y3b(1:2) = 0;
y4(1:2) = 0;
y5(1:2) = 0;

Tm(1:2) = [0,TmStp];

Lnrsd = 0;
n = 3;
while (Tm<MxTm)
    P(n) = PLwrLmt + rand*(PUpprLmt-PLwrLmt);
    [y0a, y0b, y1, y2, y3a, y3b, y4, y5] = JansenRit(y0a, y0b, y1, y2, y3a, y3b, y4, y5, C1, C2, C3, C4, a, A, b, B, e0, r, v0, P, TmStp, n, Lnrsd);
    Tm(n) = Tm(n-1) + TmStp;
    n = n + 1;
end

yr = (y1-y2);
yb = (y4-y5);

[vrbl0, NmbrPnts] = size(yr);
Tm = Tm(NmbrPnts-ceil((MxTm/2)/TmStp):NmbrPnts);
yr = yr(NmbrPnts-ceil((MxTm/2)/TmStp):NmbrPnts);
yb = yb(NmbrPnts-ceil((MxTm/2)/TmStp):NmbrPnts);
[vrbl0, NmbrPnts] = size(yr);

[f, FT] = FrrTrnsfrm(Tm, yr);

Prmtrs(1,1) = a;
Prmtrs(2,1) = A;
Prmtrs(3,1) = b;
Prmtrs(4,1) = B;
Prmtrs(5,1) = r;
Prmtrs(6,1) = e0;
Prmtrs(7,1) = v0;
Prmtrs(8,1) = C;
Prmtrs(9,1) = C2;
Prmtrs(10,1) = C3;
Prmtrs(11,1) = C4;

TF = 1/(sqrt(12))*(PUpprLmt-PLwrLmt)*TrnsfrFnctn(Prmtrs, f)/sqrt(NmbrPnts);

subplot(2,3,1);
plot(Tm, yr, 'k', 'linewidth', 2);
MnX = min(Tm);
MxX = max(Tm);
MnY = min(yr);
MxY = max(yr);
axis([MnX MxX MnY MxY]);
set(gca, 'fontsize', FntSz);
xlabel('Time (s)', 'fontsize', FntSz);
ylabel('y_1 - y_2 (mV)', 'fontsize', FntSz);

subplot(2,3,2);
title('Sigmoid', 'fontsize', FntSz);
plot(yr, yb, 'k', 'linewidth', 2);
MnX = min(yr);
MxX = max(yr);
MnY = min(yb);
MxY = max(yb);
axis([MnX MxX MnY MxY]);
set(gca, 'fontsize', FntSz);
xlabel('y_1 - y_2 (mV)', 'fontsize', FntSz);
ylabel('dy_1/dt - dy_2/dt', 'fontsize', FntSz);
title('Sigmoid', 'fontsize', FntSz);

subplot(2,3,3);
hold on;
plot(f, FT, 'k', 'linewidth', 2);
plot(f, TF, 'r', 'linewidth', 2);
legend('Fourier transform', 'Transfer function');
MnX = min(f);
MxX = max(f);
MnY = min(FT);
MxY = max(FT);
axis([MnX MxX 0 Inf]);
set(gca, 'fontsize', FntSz);
xlabel('Frequency (Hz)', 'fontsize', FntSz);
ylabel('Amplitude', 'fontsize', FntSz);


clear Tm;
clear yr;
clear yb;
clear FT;

% Linearised sigmoid function

y0a(1:2) = 0;
y0b(1:2) = 0;
y1(1:2) = 0;
y2(1:2) = 0;
y3a(1:2) = 0;
y3b(1:2) = 0;
y4(1:2) = 0;
y5(1:2) = 0;

Tm(1:2) = [0,TmStp];

Lnrsd = 1;
n = 3;
while (Tm<MxTm)
    P(n) = PLwrLmt + rand*(PUpprLmt-PLwrLmt);
    [y0a, y0b, y1, y2, y3a, y3b, y4, y5] = JansenRit(y0a, y0b, y1, y2, y3a, y3b, y4, y5, C1, C2, C3, C4, a, A, b, B, e0, r, v0, P, TmStp, n, Lnrsd);
    Tm(n) = Tm(n-1) + TmStp;
    n = n + 1;
end

yr = (y1-y2);
yb = (y4-y5);

[vrbl0, NmbrPnts] = size(yr);
Tm = Tm(NmbrPnts-ceil((MxTm/2)/TmStp):NmbrPnts);
yr = yr(NmbrPnts-ceil((MxTm/2)/TmStp):NmbrPnts);
yb = yb(NmbrPnts-ceil((MxTm/2)/TmStp):NmbrPnts);

[f, FT] = FrrTrnsfrm(Tm, yr);

subplot(2,3,4);
plot(Tm, yr, 'k', 'linewidth', 2);
MnX = min(Tm);
MxX = max(Tm);
MnY = min(yr);
MxY = max(yr);
axis([MnX MxX MnY MxY]);
set(gca, 'fontsize', FntSz);
xlabel('Time (s)', 'fontsize', FntSz);
ylabel('y_1 - y_2 (mV)', 'fontsize', FntSz);

subplot(2,3,5);
plot(yr, yb, 'k', 'linewidth', 2);
MnX = min(yr);
MxX = max(yr);
MnY = min(yb);
MxY = max(yb);
axis([MnX MxX MnY MxY]);
set(gca, 'fontsize', FntSz);
xlabel('y_1 - y_2 (mV)', 'fontsize', FntSz);
ylabel('dy_1/dt - dy_2/dt', 'fontsize', FntSz);
title('Linearised sigmoid', 'fontsize', FntSz);

subplot(2,3,6);
hold on;
plot(f, FT, 'k', 'linewidth', 2);
plot(f, TF, 'r', 'linewidth', 2);
legend('Fourier transform', 'Transfer function');
MnX = min(f);
MxX = max(f);
axis([MnX MxX 0 Inf]);
set(gca, 'fontsize', FntSz);
xlabel('Frequency (Hz)', 'fontsize', FntSz);
ylabel('Amplitude', 'fontsize', FntSz);


