%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of Examples from TOOLS FOR THE STUDY OF STABILITY AND CONVERGENCE IN SET
% DYNAMICAL SYSTEMS WITH APPLICATIONS TO FEEDBACK CONTROL
% Example: 7
% Nathalie Risso. nrisso@email.arizona.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
n = 1000;
a = 1;
X0=linspace(0.7,1.1,n)';
d1 = 0.1;
d2 = 0.05;
alpha = d2/(d1+d2);
Ball = linspace(-1,1,n);
Ball = Ball';
L=0.8;
k=-0.6;
J=10;
X=zeros(n,J);
Y=zeros(n,J);
Xe=X+max(d1,d2)*Ball;
X(:,1)=X0;
C = alpha*(1+d1*Ball)+(1-alpha)*(1+d2*Ball);
for j=1:J-1
    X(:,j+1) = a*X(:,j)+k*Xe(:,j);
    % The sensor measure is a set that contains all the possible measures
    % the sensor can actually get, at each interval.
%     z1(:,j) = X(:,j)+normrnd(0,d1,[n,1]);
%     z2(:,j) = X(:,j)+normrnd(0,d2, [n,1]);
    z1(:,j) = X(:,j).*(1+d1*Ball);
    z2(:,j) = X(:,j).*(1+d2*Ball);
    % The output is estimated using the data fusion eqn.
%    Y(:,j) =(alpha*z1(:,j)+(1-alpha)*z2(:,j));
    Y(:,j)=C.*X(:,j);
    % Dynamics of the observer:
    % Xe+ = AXe+Bu + L(Y-CXe)
    Xe(:,j+1)=a*Xe(:,j)+k*Xe(:,j)+L*((Y(:,j)-C.*(Xe(:,j))));
end
% plots
close all;
figure(1)
for i=1:J-1
 plot(0*X(:,i)+i-1,X(:,i),'b','LineWidth',3);    hold on; grid on % This is the state
 plot(0*X(:,i)+i-1,Y(:,i),'r','LineWidth',2); % This is the state estimation
 xlabel('$j$', 'Interpreter','latex');ylabel('$x_p$', 'Interpreter','latex');
end
