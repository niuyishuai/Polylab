n=5;

%% generate polynomial
disp(repmat('.',1,50));
help MPOLY_GPU
p=MPOLY_GPU(n);
p.disp;
p.sdisp;

p=MPOLY_GPU(n,2,zeros(1,n));
p.sdisp;

help MPOLY_GPU.disp
help MPOLY_GPU.sdisp

%% modify polynomial
disp(repmat('.',1,50));
p.coef = [-1; 2];
p.pow = [ 0 0 0 0 0; 1 2 3 0 1];
p.k = length(p.coef);
p.disp;
p.sdisp;

%% generate variables
disp(repmat('.',1,50));
help MPOLY_GPU.mpolyvars
x = MPOLY_GPU.mpolyvars(n);
x(1).sdisp;
x.sdisp;

%% monolist
disp(repmat('.',1,50));
help MPOLY_GPU.monolist
lst = MPOLY_GPU.monolist(n,3);
lst.sdisp;
lst.disp;

%% special polynomial matrix
%% test zeros and iszero
disp(repmat('.',1,50));
help MPOLY_GPU.zeros
help MPOLY_GPU.iszero
m = MPOLY_GPU.zeros(3);
m.sdisp;
MPOLY_GPU.iszero(m)
m = MPOLY_GPU.zeros(3,2);
m.sdisp;
MPOLY_GPU.iszero(m)
m = MPOLY_GPU.zeros(3,2,1);
m.sdisp;
MPOLY_GPU.iszero(m)
MPOLY_GPU.iszero(0)
MPOLY_GPU.iszero([0 0])
MPOLY_GPU.iszero(1)
MPOLY_GPU.iszero([1 0])

%% test ones
disp(repmat('.',1,50));
help MPOLY_GPU.ones
m = MPOLY_GPU.ones(3);
m.sdisp;
m = MPOLY_GPU.ones(3,2);
m.sdisp;
m = MPOLY_GPU.ones(3,2,1);
m.sdisp;

%% test identity
disp(repmat('.',1,50));
help MPOLY_GPU.identity
m = MPOLY_GPU.identity(3,2);
m.sdisp;

%% test diag
disp(repmat('.',1,50));
help MPOLY_GPU.diag
m = m.diag;
m.sdisp;

m = m.diag;
m.sdisp;

%% test trace
m = MPOLY_GPU.identity(3,2);
t = m.trace;
t.sdisp;

%% polynomial operations
%% test plus +
disp(repmat('.',1,50));
help MPOLY_GPU.plus
p=x(1)+x(2)+x(3);
p.sdisp;
p=p+3;
p.sdisp;
p=2+p;
p.sdisp;

%% test uplus
disp(repmat('.',1,50));
help MPOLY_GPU.uplus
q=+p;
q.sdisp;

%% test minus -
disp(repmat('.',1,50));
help MPOLY_GPU.minus
p=p-x(2);
p.sdisp;
p=2-p;
p.sdisp
p=p-2;
p.sdisp;

%% test uminus
disp(repmat('.',1,50));
help MPOLY_GPU.uminus
p=-p;
p.sdisp;

%% test simplification
disp(repmat('.',1,50));
help MPOLY_GPU.simplify
p = MPOLY_GPU(3,[1;2],[1 0 1;1 0 1]);
p.disp;
p.sdisp;
p=p.simplify;
p.disp;
p.sdisp;
% simplify with given precision
p = MPOLY_GPU(3,[1e-7;2],[1 0 1;1 0 1]);
p.sdisp;
p = p.simplify(1e-6);
p.sdisp;


%% test mtimes *
disp(repmat('.',1,50));
help MPOLY_GPU.mtimes
% polynomial scalar * polynomial scalar
p=x(1)*x(2);
p.sdisp;
% polynomial scalar * double scalar
p=p*(-2);
p.sdisp;
% double scalar * polynomial scalar
p=3*p;
p.sdisp;
% double scalar * polynomial matrix
p=2*[x(1),x(2)+x(3)];
p.sdisp;
% polynomial matrix * double scalar
p=[x(1),x(2)+x(3)]*3;
p.sdisp;
% polynomial scalar * polynomial matrix
p = x(1)*p;
p.sdisp;
% polynomial scalar * polynomial matrix
p = p*x(1);
p.sdisp;
% double matrix * polynomial matrix
p = [1;1]*p;
p.sdisp;
% polynomial matrix * polynomial matrix
p = p*p;
p.sdisp;

%% test .*
disp(repmat('.',1,50));
help MPOLY_GPU.times
p = [x(1), x(2); x(3), x(1)];
q = p.*p;
q.sdisp;

%% test mpower ^
disp(repmat('.',1,50));
help MPOLY_GPU.mpower
q=p^2;
q.sdisp;

%% test power .^
disp(repmat('.',1,50));
q=p.^2;
q.sdisp

%% mixing +, -, * and ^
disp(repmat('.',1,50));
p=5*x(1)^3 + x(1)*x(2) -2*x(3) + 3;
p.sdisp;

%% test display precision
disp(repmat('.',1,50));
p.sdisp('%.2f');
p.sdisp('%.5e');
p.sdisp;

%% polynomial matrix operations
disp(repmat('.',1,50));
m1 = [x(1)^2, x(3)^2; x(2)^2, x(4)^2];
m2 = [x(1)^2, x(3); x(2)^3, x(4)^2];
m1.disp;
m1.sdisp;
m2.disp;
m2.sdisp;

%% polynomial matrix addition m1 + m2
disp(repmat('.',1,50));
help MPOLY_GPU.plus
% polynomial matrix + polynomial matrix
m3 = m1+m2;
m3.sdisp;
% polynomial matrix + double scalar
m3 = m3 + 2;
m3.sdisp;
m3 = 1 + m3;
m3.sdisp;
% polynomial matrix  + polynomial scalar
m3 = m3 + (x(2) + x(3)^2);
m3.sdisp;
m3 = x(1) + m3;
m3.sdisp;
% polynomial matrix + double matrix;
m3 = ones(2,2) + m3;
m3.sdisp;
m3 = m3 + ones(2,2);
m3.sdisp;
% polynomial scalar + double matrix
p = x.sum;
m3 = ones(2,2);
m4 = m3+p;
m4.sdisp;
m4 = p + m4;
m4.sdisp

%% polynomial matrix substraction m1 - m2
disp(repmat('.',1,50));
help MPOLY_GPU.minus
m3 = m3 - m2;
m3.sdisp;

%% unary minus -m3
disp(repmat('.',1,50));
help MPOLY_GPU.uminus
m4 = -m3;
m4.sdisp;

%% elementwise times .*
disp(repmat('.',1,50));
help MPOLY_GPU.times
% matrix .* matrix
m5 = m3.*m2;
m3.sdisp;
m2.sdisp;
m5.sdisp;

% matrix .* scalar
m5 =m1.*3;
m5.sdisp;

%% matrix mtimes *
%% double scalar * polynomial matrix
disp(repmat('.',1,50));
help MPOLY_GPU.mtimes
m1 = [x(1),x(2);x(3),x(1)];
m6 = 3*m1;
m7 = m1*3;
m1.sdisp;
m6.sdisp;
m7.sdisp;

%% polynomial scalar * polynomial matrix
m1 = [x(1),x(2);x(3),x(1)];
m6 = x(1)*m1;
m7 = m1*x(1);
m1.sdisp;
m6.sdisp;
m7.sdisp;

%% polynomial matrix * polynomial matrix
m8 = m1*m2;
m1.sdisp;
m2.sdisp;
m8.sdisp;

%% polynomial matrix * double vector
m8 = m8*ones(2,1);
m8.sdisp;
m8 = ones(1,2)*m8;
m8.sdisp;

%% test elementwise power .^
disp(repmat('.',1,50));
help MPOLY_GPU.power
m9 = m1.^3;
m9.sdisp;

%% test mpower ^
disp(repmat('.',1,50));
help MPOLY_GPU.mpower
m9 = m1^3;
m9.sdisp;

%% polynomial features
%% get degree
disp(repmat('.',1,50));
help MPOLY_GPU.degree
m9.degree

%% get monomials
disp(repmat('.',1,50));
help MPOLY_GPU.mono
lst=m9(1).mono();
lst.sdisp
lst=m9(1).mono(1)
lst.sdisp

%% get coefficients
disp(repmat('.',1,50));
help MPOLY_GPU.coefficients
coef = m9(1).coefficients;
disp(coef);
[~,monolst] = m9(1).coefficients;
monolst.sdisp;

%% transpose and ctranspose
disp(repmat('.',1,50));
m10 = m1';
m1.sdisp;
m10.sdisp;

%% == and ~=
% eq
disp(repmat('.',1,50));
help MPOLY_GPU.eq
help MPOLY_GPU.ne
m11=m10;
m11==m10

m11=m10+1;
m11==m10

% ne
m11=m10;
m11~=m10

m11=m10+1;
m11~=m10

%% calcul diff
disp(repmat('.',1,50));
help MPOLY_GPU.jacobian
p = x(1)+2*x(2) - x(1)*x(2) + 3*x(1)^2 - x(2)^2;
grad = p.jacobian;
p.sdisp;
grad.sdisp;

hess = grad.jacobian;
hess.sdisp;

%% evaluate polynomial
x=ones(n,1);
fval = hess.eval(x);
disp(fval);

%% test mpoly2yalmip
help mpoly2yalmip
x_yalmip = sdpvar(n,1);
m10_yalmip=mpoly2yalmip(m10,x_yalmip);
sdisplay(m10_yalmip)

%% test yalmip2mpolygpu
help yalmip2mpolygpu
m10_mpoly=yalmip2mpolygpu(m10_yalmip,x_yalmip);
m10_mpoly.sdisp;

%% test cpulize
help MPOLY_GPU.cpulize
m10_cpu=MPOLY_GPU.cpulize(m10_mpoly);
struct(m10_cpu)

%% test gpulize
help MPOLY_GPU.gpulize
m10_gpu=MPOLY_GPU.gpulize(m10_cpu);
struct(m10_gpu)

%% compare gpu and cpu MPOLY
n=5;
d=3;
m=MPOLY.monolist(n,d);
p=sum(m);

tic
p1=simplify(p^4);
t1=toc;
fprintf('CPU version: %f seconds\n',t1);

pg=MPOLY_GPU.gpulize(p);

tic
p2=simplify(pg^4);
t2=toc;
fprintf('GPU version: %f seconds\n',t2);


fprintf('Check results between GPU and CPU\n');
sdisp(simplify(MPOLY_GPU.cpulize(p2)-p1));