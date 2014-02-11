clear all
close all
clc

% run /home/mzanon/mario_data/Optimization_software/cvx/cvx_setup

syms xc_0 x_0 y_0 dxc_0 dx_0 dy_0 real
syms xc_1 x_1 y_1 dxc_1 dx_1 dy_1 real
syms lambda_0 lambda_1 real
syms u_0 u_1 real
% syms dt real

dt = 0.1;
N = 10;

s_0 = [xc_0 x_0 y_0 dxc_0 dx_0 dy_0].';
s_1 = [xc_1 x_1 y_1 dxc_1 dx_1 dy_1].';

V_0 = [s_0; lambda_0; u_0];
V_1 = [s_1; lambda_1; u_1];

M = 1;
m = 0.1;
g = 9.81;
l = 0.5;


H = [diag([M m m]), [xc_1-x_1; x_1-xc_1; y_1];
     [xc_1-x_1; x_1-xc_1; y_1].', 0];

rhs = [0; u_0; -m*g; -y_1^2 - (xc_1-x_1)^2];

f = [(xc_1-xc_0)/dt - xc_1;
     (x_1-x_0)/dt - x_1;
     (y_1-y_0)/dt - y_1;
     M*(dxc_1-dxc_0)/dt + (xc_1-x_1)*lambda_1 - u_0;
     m*(dx_1-dx_0)/dt  + (x_1-xc_1)*lambda_1;
     m*(dy_1-dy_0)/dt  +      y_1*lambda_1 + m*g;
     (xc_1-x_1)*(dxc_1-dxc_0)/dt + (x_1-xc_1)*(dx_1-dx_0)/dt + y_1*(dy_1-dy_0)/dt + dy_1^2 + (dxc_1-dx_1)^2;];


for k = 1:length(f)
    F(k) = subs(f(k),[V_0;V_1],0*[V_0;V_1]);
    q(:,k) = subs(jacobian(f(k),[V_0;V_1]).',[V_0;V_1],0*[V_0;V_1]);
    Q(:,:,k) = jacobian(jacobian(f(k),[V_0;V_1]),[V_0;V_1]);
    Qconst(:,:,k) = [2*F(k), q(:,k).';
                     q(:,k), Q(:,:,k)];
end

V = [1; V_0; V_1]*[1; V_0; V_1].';
for k = 1:length(f)
    simplify(f(k) - trace(Qconst(:,:,k)*V)/2)
end

Xsize = 1 + length(V_0)*2;

Qunitary = zeros(Xsize,Xsize);
Qunitary(1,1) = 1;


theta = 0;
xc = 0;
x =  xc + l*sin(theta);
y =      -l*cos(theta);

Xinit = [xc x y 0 0 0];
% Vinit = [Xinit 0 0];

Q_init = [diag([ones(1,6) 0 0]) zeros(8)];

Xref = [0 0 -l 0 0 0].';
Uref = 0;
ref = [1; Xref; 0; Uref; Xref; 0; Uref];


Wx = ([1 10 10 1e-2 1e-2 1e-2]);
Wu = ([1e-2]);
Wv = diag([Wx 0 Wu]);
Qcost = blkdiag(0,Wv,0*Wv);
Qcost(1,2:1+length(V_0)) = Wv*[Xref; 0; Uref;];
Qcost(2:1+length(V_0),1) = (Wv*[Xref; 0; Uref;]).';

WvT = diag([Wx 1 Wu]);
QcostT = blkdiag(0,0*Wv,WvT);

a = 1e0;


cvx_begin
    
    variable X(Xsize,Xsize,N)
    
    cost = 0;
    for k = 1:N
        cost = cost + trace(Qcost*X(:,:,k)) + a*norm_nuc(X(:,:,k)/length(X(:,:,k))^2);

        X(:,:,k) == semidefinite(Xsize);
    end
    % Weight on the terminal state
    cost = cost + trace(QcostT*X(:,:,end));
    
    minimize( cost );
    subject to
    % System dynamics at each time instant ...
    for k = 1:N
        % For every state
        for j = 1:length(s_0)
            trace(double(Qconst(:,:,j))*X(:,:,k) ) == 0;
        end
    end
    
    % Connect the variables at every node ...
    for k = 2:N
        X(1,2:length(V_0)+1,k) == X(1,2+length(V_0):end,k-1);
    end

    % Initial condition
    for j = 1:length(s_0)
        X(1,1+j,1) - Xinit(j) == 0;
    end
    % Final condition
    
    for k = 1:N
        X(1,1,k) == 1;
%         trace(Qunitary*X(:,:,k)) == 1;
    end
    
cvx_end

check = [];
for k = 1:N
    check = blkdiag(check,X(:,:,k));
end

s = svd(check);
s(1)/s(2)


state = [];
lambda = [];
control = [];
v = [];
for k = 1:N
    v = [v; X(1,2:2+length(V_0),k)];
    state = [state; X(1,2:1+length(s_0),k)];
    lambda = [lambda; X(1,2+length(s_0),k)];
    control = [control; X(1,1+length(V_0),k)];
end

v
