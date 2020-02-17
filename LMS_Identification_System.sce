 //identificação de sistemas utilizando LMS
clear;
clc;

rand("normal");
N = 2000;
hn = [1.2 0.8 0.6];
sigma1 = 0.2;

E = rand(1,N);
xn = filter(hn,1,E)';
W = [1 2 3];
n = sigma1*rand(1,N);

rx = xcorr(xn,2,'biased');
Rx = [rx(3:5) rx(2:4) rx(1:3)];
tr = Rx(1,1)+Rx(2,2)+Rx(3,3);
[Q,lambda]=spec(Rx);
lambdamax=max(lambda);
// Computacao iterativa
c = 0.1
mi = c/(2*lambdamax+tr);
wn = [0 0 0]'; // filtro unitario inicial
d = filter(W,1,xn)+n';
w = zeros(3,N);

xb = [0 0 0]';
for n=1:length(xn)
    xb = [xn(n);xb(1:2)];
    yn = wn'*xb;
    y(n) = yn;
    e(n) = d(n)-y(n);
    // calculo do grad estimado
    grade = -2*e(n)*xb;
    wn = wn - mi*grade;
    w(:,n) = wn
end

plot2d(w');
