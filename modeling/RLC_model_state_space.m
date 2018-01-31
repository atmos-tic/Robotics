R = 1;
L = 1;
C = 1;

A = [0, 1; 1/(L*C), R/L]
b = [0, 1/(L*C)]'
c = [1, 0]
d = 0

S = ss(A,b,c,d);
impulseplot(S(1));