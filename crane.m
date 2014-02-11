clear all
close all
clc

% run /home/mzanon/mario_data/Optimization_software/cvx/cvx_setup

syms xc_0 x_0 y_0 dxc_0 dx_0 dy_0 real
syms xc_1 x_1 y_1 dxc_1 dx_1 dy_1 real
syms lambda_0 lambda_1 real
syms u_0 u_1 real
% syms dt real

dt = 0.01;
N = 100;

s_0 = [xc_0 x_0 y_0 dxc_0 dx_0 dy_0].';
s_1 = [xc_1 x_1 y_1 dxc_1 dx_1 dy_1].';

V_0 = [s_0; lambda_0; u_0];
V_1 = [s_1; lambda_1; u_1];

M = 1;
m = 0.1;
g = 9.81;
l = 0.5;


% H = [diag([M m m]), [xc_1-x_1; x_1-xc_1; y_1];
%      [xc_1-x_1; x_1-xc_1; y_1].', 0];
% 
% H = blkdiag(eye(3),H);
% 
% DV = (V_1 - V_0)/dt;
% DV = [DV(1:6);DV(7)*dt];
% 
% rhs = [0; 0; 0; u_0; 0; -m*g; -y_1^2 - (xc_1-x_1)^2];
% 
% H*DV - rhs;

scale = 1e2;

cstr = 0.5*( (xc_1-x_1)^2 + y_1^2 - l^2 );
dcstr = (dxc_1-dx_1)*(xc_1-x_1) + dy_1*y_1;

p = 0;

f = [(xc_1-xc_0)/dt - dxc_1;
     (x_1-x_0)/dt - dx_1;
     (y_1-y_0)/dt - dy_1;
     M*(dxc_1-dxc_0)/dt + (xc_1-x_1)*scale*lambda_1 - u_0;
     m*(dx_1-dx_0)/dt   + (x_1-xc_1)*scale*lambda_1;
     m*(dy_1-dy_0)/dt   +        y_1*scale*lambda_1 + m*g;
     (xc_1-x_1)*(dxc_1-dxc_0)/dt + (x_1-xc_1)*(dx_1-dx_0)/dt + y_1*(dy_1-dy_0)/dt + dy_1^2 + (dxc_1-dx_1)^2  + p^2*cstr + p*dcstr;];


% f = [(xc_1-xc_0)/dt;
%      (x_1-x_0)/dt;
%      (y_1-y_0)/dt - dy_1;
%      M*(dxc_1-dxc_0)/dt;
%      m*(dx_1-dx_0)/dt;
%      m*(dy_1-dy_0)/dt  +      y_1*scale*lambda_1 + m*g;
%      y_1*(dy_1-dy_0)/dt + dy_1^2;];
% f(end) = y_1*(dy_1-dy_0)/dt + dy_1^2;


for k = 1:length(f)
    F(k) = subs(f(k),[V_0;V_1],0*[V_0;V_1]);
    q(:,k) = subs(jacobian(f(k),[V_0;V_1]).',[V_0;V_1],0*[V_0;V_1]);
    Q(:,:,k) = jacobian(jacobian(f(k),[V_0;V_1]),[V_0;V_1]);
    Qconst(:,:,k) = [2*F(k), q(:,k).';
                     q(:,k), Q(:,:,k)];
end

V = [1; V_0; V_1]*[1; V_0; V_1].';
% for k = 1:length(f)
%     simplify(f(k) - trace(Qconst(:,:,k)*V)/2)
% end

Xsize = 1 + length(V_0)*2;

Qunitary = zeros(Xsize,Xsize);
Qunitary(1,1) = 1;


theta = 0;
xc = 0.2;
x =  xc + l*sin(theta);
y =      -l*cos(theta);

Xinit = [xc x y 0 0 0];
% Vinit = [Xinit 0 0];

Q_init = [diag([ones(1,6) 0 0]) zeros(8)];

Xref = [0.2 0.2 l 0 0 0].';
Uref = 0;
ref = [Xref; 0; Uref;];


Wx = ([1 10 10 1e-2 1e-2 1e-2]);
Wu = ([1e-2]);
Wv = diag([Wx 0 Wu]);

c = (V_0-ref).'*Wv*(V_0-ref);
C = subs(c,[V_0;V_1],0*[V_0;V_1]);
qc = subs(jacobian(c,[V_0;V_1]).',[V_0;V_1],0*[V_0;V_1]);
Qc = jacobian(jacobian(c,[V_0;V_1]),[V_0;V_1]);
Qcost = double([2*C, qc.';
          qc, Qc]/2);

% Qcost = blkdiag(0,Wv,0*Wv);
% Qcost(1,2:1+length(V_0)) = Wv*[Xref; 0; Uref;];
% Qcost(2:1+length(V_0),1) = (Wv*[Xref; 0; Uref;]).';

WvT = diag([Wx 1 Wu]);
QcostT = blkdiag(0,0*Wv,WvT);

a = 1e4;

disp('Prepared the system')

cvx_begin

    cvx_precision best

    variable X(Xsize,Xsize,N)
    
    cost = 0;
    for k = 1:N
        cost = cost + trace(Qcost*X(:,:,k));

        X(:,:,k) == semidefinite(Xsize);
    end
    % Weight on the terminal state
    cost = cost + trace(QcostT*X(:,:,end));
    cost2 = cost;
    for k = 1:N
        cost = cost + a*norm_nuc(X(:,:,k))/length(X(:,:,k))^2;
    end
    
    disp('Created cost function')
    
%     string = 'XX=blkdiag(';
%     for k = 1:N
%         string = [string,'X(:,:,',num2str(k),'),'];
%     end
%     string = [string(1:end-1),');'];
%     eval(string)
% 
%     cost = cost + a*norm_nuc(XX)/length(XX)^2;
    
    minimize( cost );
    subject to
    % System dynamics at each time instant ...
    for k = 1:N
        % For every state
        for j = 1:length(f)
            trace(double(Qconst(:,:,j))*X(:,:,k) ) == 0;
        end
    end
    
    disp('Created dynamic constraints')
    
    % Connect the variables at every node ...
    for k = 2:N
        X(1,2:length(V_0)+1,k) == X(1,2+length(V_0):end,k-1);
    end
    
    disp('Created connecting constraints between stages')

    % Initial condition
    for j = 1:length(s_0)
        X(1,1+j,1) - Xinit(j) == 0;
    end
    % Final condition
    for j = 1:length(s_0)
        X(1,1+j,end) - Xref(j) == 0;
    end
    
    disp('Created initial and final constraints')
    
    for k = 1:N
        X(1,1,k) == 1;
%         trace(Qunitary*X(:,:,k)) == 1;
    end
    disp('Created 1 constraints for each stage')
    
cvx_end

% s = svds(XX,length(XX));
% s(1)/s(2)
for k = 1:N
    s = svd(X(:,:,k));
    s(1)/s(2)
end

state = [];
lambda = [];
lambda1 = [];
control = [];
v = [];
v1 = [];
for k = 1:N
    v = [v; X(1,2:1+length(V_0),k)];
    v1 = [v1; X(1,2+length(V_0):end,k)];
    state = [state; X(1,2:1+length(s_0),k)];
    lambda = [lambda; X(1,2+length(s_0),k)];
    lambda1 = [lambda1; X(1,end-1,k)];
    control = [control; X(1,1+length(V_0),k)];
end


v
v1
lambda1
cost2


z = [Xref;m*g/0.5/scale;0];
z = [1;z;z];
Z = z*z.';

% check = [];
% Check = [];
% CH = [];
% for k = 1:N
%     check0 = [];
%     Check0 = [];
%     CH0 = [];
%     for j = 1:length(f)
%         check0 = [check0, subs(f(j),[V_0; V_1],X(2:end,1,k))];
%         Check0 = [Check0, subs(trace(double(Qconst(:,:,j))*V )/2,[V_0; V_1],X(2:end,1,k))];
%         CH0 = [CH0, trace(double(Qconst(:,:,j))*X(:,:,k) )/2];
%     end
%     check = [check;check0];
%     Check = [Check;Check0];
%     CH = [CH;CH0];
% end
% check
% Check
% CH


%% Plot

figure(1)
clf
plot(v(:,1),v(:,1)*0,'.k')
hold on
plot(v(:,2),v(:,3),'.b')

figure(2)
clf
subplot(3,1,1)
hold on
plot((0:N-1)*dt,v(:,1),'k')
plot((0:N-1)*dt,v(:,2),'b')
subplot(3,1,2)
plot((0:N-1)*dt,v(:,3),'b')

subplot(3,1,3)
stairs((0:N-1)*dt,v(:,end))

