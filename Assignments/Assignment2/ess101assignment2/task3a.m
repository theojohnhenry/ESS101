clc,clear, clf
rng(2)

syms x real

N = 100;
anoll = -2.923;
bnoll = 7.18; 

ilin = linspace(1,N,N)
xlin = unifrnd(0,50,[N, 1])
elin = normrnd(0,3.8, [N,1])

ylin = anoll + bnoll*xlin+elin;

scatter(xlin,ylin)
hold on 

bigphi = [ones(N,1), xlin];
thetahat = inv(bigphi'*bigphi)*bigphi'*ylin

xplutt = linspace(1, 50, N); 
yplutt = thetahat(1) + thetahat(2)*xplutt;     

plot(xplutt,yplutt,'r','LineWidth',1)

yhat = bigphi * thetahat;
epsilon = mean((ylin - yhat).^2)