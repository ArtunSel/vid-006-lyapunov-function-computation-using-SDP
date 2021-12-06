clear all,close all,clc;
%%
% this is for SDPT3
current_dir=pwd;
cd('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\SDPT3-4.0');
run('Installmex.m')
run('startup.m')
cd(current_dir);
% cd('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\SDPT3-4.0\Solver')
% run('sqlpdemo.m')

% this is for YALMIP
addpath(genpath('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\YALMIP-master'))
yalmip('clear');

% this is for SOSTOOLS
cd('C:\Users\Artun\Documents\MATLAB\SOSTOOLS.303');
run('addsostools.m')

% go back to the original directory
cd(current_dir);
clear all
%% yalmip 6-hump camel optimization proble
clear all,close all,clc;yalmip('clear');
sdpvar x1 x2
sdpvar gamma
p=4*x1^2-2.1*x1^4+(1/3)*x1^6+x1*x2-4*x2^2+4*x2^4;

F=[sos(p-gamma)];



ops = sdpsettings('solver','sdpt3');
ops.dualize=[0];
ops.sos.newton=[1];
ops.sos.congruence=[1];
ops.sos.numblkdg=[1e-6];
ops.sdpt3.maxit=100;
sol = optimize(F,[-gamma],ops);
sol.info

value(gamma)

%% matlab plot of the six-hump camel function
aa1=2;
aa2=1;
% aa1=100;
% aa2=100;
x1=linspace(-aa1,aa1,50);
x2=linspace(-aa2,aa2,50);

[X1,X2] = meshgrid(x1,x2);
% Z = 4*X1.^2+X1.*X2-4*X2.^2-2.4*X1.^4+4*X2.^4+(1/3)*X1.^6;
% Z = 4*X1.^2-2.1*X1.^4+(1/3)*X1.^6+X1.*X2-4*X2.^2+4*X2.^4;
Z = 4*X1.^2-2.1*X1.^4+(1/3)*X1.^6+X1.*X2-4*X2.^2+4*X2.^4;
surf(X1,X2,Z)
xlabel('X1');
ylabel('X2');
%%



%% YALMIP sos decomposition of a p(x,y) vid-1
clear all,close all,clc;yalmip('clear');
x = sdpvar(1,1);y = sdpvar(1,1);
% p = (1+x)^4 + (1-y)^2;
p = 2*x^4+5*y^4-x^2*y^2+2*x^3*y;

degree(p)
ceil(degree(p)/2)

F=[sos(p)];

ops = sdpsettings('solver','sdpt3');
ops.dualize=[0];
ops.sos.newton=[1];
ops.sos.congruence=[1];
ops.sos.numblkdg=[1e-6];
ops.sdpt3.maxit=100;

[sol,m,QX,residuals,everything] = solvesos(F,[],ops,[],monolist([x y],2,1));
sol.info
checkset(F)
sdisplay(m{1})
min(eig(QX{1}))



[c,v] = coefficients(m{1}'*QX{1}*m{1},[x y]);        sdisplay([c,v])
%% ----------------------


%----------------------
h = sosd(F);
sdisplay(h)

clean(p-h'*h,1e-6)

chol(QX{1})




%% YALMIP sos decomposition of a polynomail-matrix P(x,y) vid-2 "READY"
clear all,close all,clc;yalmip('clear');
x = sdpvar(1,1);y = sdpvar(1,1);
% p = 4*x^2-2.1*x^4+(1/3)*x^6+x*y-4*y^2+4*y^4;
P = [1+x^2 -x+y+x^2;-x+y+x^2 2*x^2-2*x*y+y^2];
degree(P)



ops = sdpsettings('solver','sdpt3');
ops.dualize=[0];
ops.sos.newton=[1];
ops.sos.congruence=[1];
ops.sos.numblkdg=[1e-6];
ops.sdpt3.maxit=100;

F=[sos(P)];
[sol,m,QX,residuals,everything] = solvesos(F);
sol.info

sdisplay(m{1})
min(eig(QX{1}))

sdisplay(m{1}'*QX{1}*m{1})
temp1=m{1}'*QX{1}*m{1};
temp1=clean(temp1,1e-5)

sdisplay(temp1) % which is equalt to the poly-matrix P(x,y)










%% YOUTUBE VIDEO problem YALMIP
clear all,close all,clc
a1 = sdpvar(1,1);a2 = sdpvar(1,1);
eps1=1e-5;
F=[diag([a1;a2])>=eps1*eye(2)];
F=[F;[-2*a2,0,2*a1-a2;0,-4*a2,0;2*a1-a2,0,-2*a1]<=eps1*eye(3)];
F=[F;1 <= [a1 a2] <= 5];
ops = sdpsettings('solver','sdpt3');
sol=optimize(F,[],ops);
sol.info

a10=value(a1);
a20=value(a2);


min(eig(diag([a10;a20])))
eig([-2*a20,0,2*a10-a20;0,-4*a20,0;2*a10-a20,0,-2*a10])
%% Vdot computation
syms x1 x2 real
syms a1 a2 real

f1=-x1+2*x2^2;
f2=-x2+-x1*x2-2*x2^3;
V=a1*x1^2+a2*x2^2;
Vdot=jacobian(V,[x1,x2])*[f1;f2]
Vdot=expand(Vdot);

[C1,T1] = coeffs(Vdot,[x1,x2])
max_deg=max(polynomialDegree(T1))






%