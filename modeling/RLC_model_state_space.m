R = 1;
L = 1;
C = 1;

A = [0, 1; -1/(L*C), -R/L];
b = [0; 1/(L*C)];
Bu = [1;1];
c = [1, 0];
d = 0;

Q=1; R=10;   %雑音
N = 50;

v = randn(N,2)*sqrtm(Q);    %システム雑音
w = randn(N,1)*sqrtm(R);    %観測雑音
x = zeros(N,2); y = zeros(N,1); %記憶領域の確保
u = ones(N,2);

y(1) = c*x(1,:)'+w(1);
 for k = 2:N
        x(k,:) = A*x(k-1,:)' + b*u(k-1) + Bu*v(k-1);
        y(k) = c*x(k,:)'+w(k);
 end
 
xhat = zeros(N,2);
P = 0; xhat(1,:) = 0;
for k = 2:N
    [xhat(k,:),P,G] = kf(A,b,0,c,Q,R,0,y(k),xhat(k-1,:),P);
end
figure(1),clf
plot(1:N,y,'k:',1:N,x(:,1),'r--',1:N,xhat(:,1),'b-')
xlabel('No. of samples')
legend('measured','true','estimate')

%S = ss(A,b,c,d);
%G = tf(S);
%impulseplot(G);

%%線形カルマンフィルタのfunction文
function [xhat_new,P_new, G] = kf(A,B,Bu,C,Q,R,u,y,xhat,P)
% KF 線形カルマンフィルタの更新
% [xhat_new,P_new, G] = kf(A,B,Bu,C,Q,R,u,y,xhat,P)
% 線形カルマンフィルタの推定値更新を行う
% 引数:
%   A,B,h,C: 対象システム
%           x(k+1) = Ax(k) + Bv(k) + BU u(k)
%           y(k) = C'x(k) + w(k)
%           のシステム行列
%   Q,R: 雑音v,wの共分散行列. v,w は正規性白色雑音で
%           E[v(k)] = E[w(k)] = 0
%           E[v(k)'v(k)] = Q, E[w(k)'w(k)] = R
%       であることを想定
%   u: 状態更新前時点での制御入力 u(k-1)
%   y: 状態更新後時点での観測出力 y(k)
% 戻り値:
%   xhat_new: 更新後の状態推定値 xhat(k)
%   P_new: 更新後の誤差共分散行列 P(k)
%   G: カルマンゲイン
% 参考:
%   非線形システムへの拡張: EKF, UKF
% 列ベクトルに整形
 xhat = xhat(:); u = u(:); y = y(:);
%事前推定値
 xhatm = A*xhat + Bu*u; %状態
 Pm = A*P*A' + B*Q*B';   %誤差共分散
%カルマンゲイン行列
 G = Pm*C'/(C*Pm*C'+R);
% 事後推定値
 xhat_new = xhatm + G * (y-C*xhat);    %状態
 P_new = (eye(size(A))-G*C) * Pm;  %誤差共分散
end
 