R = 1;
L = 1;
C = 1;

A = [0, 1; -1/(L*C), -R/L];
b = [0; 1/(L*C)];
Bu = [1;1];
c = [1, 0];
d = 0;

Q=1; R=10;   %�G��
N = 10;

v = randn(N,2)*sqrtm(Q);    %�V�X�e���G��
w = randn(N,1)*sqrtm(R);    %�ϑ��G��
x = zeros(N,2); y = zeros(N,1); %�L���̈�̊m��
u = ones(N,2);

y(1) = c*x(1,:)';%+w(1);
 for k = 2:N
        x(k,:) = A*x(k-1,:)' + b*u(k-1) + Bu*v(k-1);
        y(k) = c*x(k,:)';%+w(k);
 end
 
xhat = zeros(N,2);
P = 0; xhat(1,:) = 0;
for k = 2:N
    [xhat(k,:),P,G] = kf(A,b,Bu,c,Q,R,u(k),y(k),xhat(k-1,:),P);
end
figure(1),clf
plot(1:N,y,'k:',1:N,x(:,1),'r--',1:N,xhat(:,1),'b-')
xlabel('No. of samples')
legend('measured','true','estimate')

%S = ss(A,b,c,d);
%G = tf(S);
%impulseplot(G);

%%���`�J���}���t�B���^��function��
function [xhat_new,P_new, G] = kf(A,B,Bu,C,Q,R,u,y,xhat,P)
% KF ���`�J���}���t�B���^�̍X�V
% [xhat_new,P_new, G] = kf(A,B,Bu,C,Q,R,u,y,xhat,P)
% ���`�J���}���t�B���^�̐���l�X�V���s��
% ����:
%   A,B,h,C: �ΏۃV�X�e��
%           x(k+1) = Ax(k) + Bv(k) + BU u(k)
%           y(k) = C'x(k) + w(k)
%           �̃V�X�e���s��
%   Q,R: �G��v,w�̋����U�s��. v,w �͐��K�����F�G����
%           E[v(k)] = E[w(k)] = 0
%           E[v(k)'v(k)] = Q, E[w(k)'w(k)] = R
%       �ł��邱�Ƃ�z��
%   u: ��ԍX�V�O���_�ł̐������ u(k-1)
%   y: ��ԍX�V�㎞�_�ł̊ϑ��o�� y(k)
% �߂�l:
%   xhat_new: �X�V��̏�Ԑ���l xhat(k)
%   P_new: �X�V��̌덷�����U�s�� P(k)
%   G: �J���}���Q�C��

% �Q�l:
%   ����`�V�X�e���ւ̊g��: EKF, UKF
% ��x�N�g���ɐ��`
 xhat = xhat(:); u = u(:); y = y(:);
%���O����l
 xhatm = A*xhat + B*u; %���
 Pm = A*P*A' + Bu*Q*Bu';   %�덷�����U
%�J���}���Q�C���s��
 G = Pm*C'/(C*Pm*C'+R);
% ���㐄��l
 xhat_new = xhatm + G * (y-C*xhat);    %���
 P_new = (eye(size(A))-G*C) * Pm;  %�덷�����U
end
 