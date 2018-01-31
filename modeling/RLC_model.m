R = 1;
L = 1;
C = 1;

G = tf(1,[C*L, C*R, 1]);
S = ss(G)
stepplot(S);
