clc,clear, clf
rng(2)

syms x real

N = 100;
anoll = -2.923;
bnoll = 7.18; 
cnoll = 2.8;

ilin = linspace(1,N,N)
xlin = unifrnd(0,50,[N, 1])
xsquaredlin = xlin.^2
elin = normrnd(0,12.8, [N,1])

ylin = anoll + bnoll*xlin+elin + cnoll*xsquaredlin

scatter(xlin,ylin)
hold on 

bigphi = [ones(N,1), xlin, xsquaredlin]
thetahat = inv(bigphi'*bigphi)*bigphi'*ylin

xplutt = linspace(1, 50, N); 
yplutt = thetahat(1) + thetahat(2)*xplutt + thetahat(3)*xplutt.^2;        

plot(xplutt,yplutt,'r','LineWidth',1)

yhat = bigphi * thetahat;
epsilon = mean((ylin - yhat).^2)