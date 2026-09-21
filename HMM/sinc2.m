function y = sinc2(x)

y = sin(pi*x)./(pi*x);
y(x == 0) = 1;