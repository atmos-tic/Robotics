R = 1;
L = 1;
C = 1;

G = tf(1,[C*L, R*L, 1]);
stepplot(G);
impulseplot(G);