clear
close
t = (-5:0.02:5);
B = [1, 0.5, 0.3, 0.2, 0.1]';
T = 1;
ccdf = @(x) 1-cdf(makedist('Normal'), x);
g = 1/(2*T)*(ccdf(2*pi*B*(t-T/2)/sqrt(log(2))) ...
            -ccdf(2*pi*B*(t+T/2)/sqrt(log(2))));

plot(t,g)
ylim([0, 0.55])
grid on
legend('BT=1', 'BT=0.5', 'BT=0.3', 'BT=0.2', 'BT=0.1')