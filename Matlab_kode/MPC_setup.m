Model = 2;

if (Model == 1)
    % Discrete double integrator
    A = [1 1;0 1];
    B = [1;0.3];  %'Double integrator from C&CE paper with Florin
    C = eye(2);

    nx = 2;
    ny = 2;
    nu = 1;
    n = 10;

    sys = ss(A,B,C,0,1);

    qy = eye(2); %Weight on y
    Q = C'*qy*C
    R = 1; %Weight on u

    umin = -1;
    umax = 1;
    ymin = -5*[1;1];   %Constraints on _states_ since C =eye(nx)
    ymax = 5*[1;1];

end

if (Model == 2)
    % 4-state example from paper with Florin, example 2
    A = [0.928 0.002 -0.003 -0.004; 0.041 0.954 0.012 0.006;-0.052 -0.046 0.893 -0.003; -0.069 0.051 0.032 0.935];
    B = [0 0.336;0.183 0.007;0.090 -0.009;0.042 0.012];  
    C = [0 0 -0.098 0.269;0 0 0.080 0.327];

    nx = 4;
    ny = 2;
    nu = 2;
    n = 30;

    sys = ss(A,B,C,0,1); %creating state space model

    qy = eye(2); %Weight on y
    Q = C'*qy*C;
 %   Q = Q + daug(0.01*eye(2),zeros(2,2));
    R = eye(nu); %Weight on u
 %   R = 5*eye(nu);

    umin = [-1;-1];
    umax = [1;1];
    ymin = [-1;-1];   
    ymax = [1;1];
end

if (Model == 3)
    %Skogestad 82-state distillations Column A.  The model is scaled, to
    %maximum 1 in inputs and outputs(offsets).
    %Discretized eith 0.1 min sample time.
    %Outputs: Top composition, Bottom composition, Top accumulator holdup,
    %   Bottom holdup
    %Inputs: Top reflux, Bottom boilup, Top distillate rate, Bottom product
    %rate
    %Normally we would close the two level loops and have only a 2x2
    %system, but here we will try controlling it directly with MPC as a 4x4
    %system.  Hence also the somewhat short sample period for composition
    %control (needed for level control).
    load colA.mat
    
    %Unable to calculate terminal set for 82-state system, so do model
    %reduction to 25 states
    [Gb,sg] = balreal(Gcont);
    elim = (sg<1e-5);    %1e-5
    Gr = modred(Gb,elim);
    
    Grd = c2d(Gr,0.1);
    
 %   Grd = c2d(Gcont,0.1);
    
    [A,B,C,D] = ssdata(Grd);
    umin = -ones(4,1);
    umax = ones(4,1);
    ymax = ones(4,1);
    ymin = -ones(4,1);
    
    %nx = 82;
    nx = 25;
    nu = 4;
    ny = 4;
    n = 30;
    Q = C'*eye(ny)*C;
    R = eye(nu);
    sys = Grd;
end

%K is the gain matrix 
[K,S,E]=lqr(sys,Q,R);

%Hva skjer her?
%kron returns the Kronecker tensor product
L=kron(eye(n-1),Q);
Qh = daug(kron(eye(n-1),Q),S);
Rh = kron(eye(n),R); 



%variable orgainization:[x(1:np);u(0;np-1)]

%ymax constraints
Hy = kron(eye(n-1),C);
hy = kron(ones(n-1,1),ymax); %making a ymax vector (repeating ymax n-1 times) 
%ymin constraint
F = kron(eye(n-1),-C);
Hy = [Hy;kron(eye(n-1),-C)];
hy = [hy;kron(ones(n-1,1),-ymin)]; %adding the ymin constraints to hy (repeating ymin n-1 times) 
%Terminal constraints
tic
[Hn,hn] = termset(A,B,C,K,ymax,ymin,umax,umin);

if Model == 1
    save -ascii 'A.txt' Hn;
    save -ascii 'B.txt' hn;
end
if Model == 2
    save -ascii 'A2.txt' Hn;
    save -ascii 'B2.txt' hn;
end
toc
Hy = daug(Hy,Hn); %dowload robust control toolbox
hy = [hy;hn];
%umax constraints
Hu = eye(n*nu);
hu = kron(ones(n,1),umax);
%umin constraints 
Hu = [Hu;-eye(n*nu)];
hu = [hu;kron(ones(n,1),-umin)];


%above: specified the model and the constraints, now comes the solution to
%the problem


v = ones(n-1,1); %generating vector [1; 1; ...] with n-1 cols
Ia = eye(n*nx)-kron(diag(v,-1),A); %Identitiy matrix - [constucted matrix containing A]
Iai = inv(Ia);
Bh = kron(eye(n),B);
%Ch = kron(eye(n),C);

O =[1;zeros(n-1,1)];
A0 = kron([1;zeros(n-1,1)],A);

% xvec = inv(Ia)*A0*x0+inv(Ia)*Bh*uvec
% Gu \leq Ss*x_0+W
Gy = Hy*Iai*Bh;
Ssy = -Hy*Iai*A0;
Wy = hy;
                %Not very efficient this ;-)
Gu = Hu;
Ssu = zeros(size(hu,1),nx);
Wu = hu;

G = [Gy;Gu];
Ss = [Ssy;Ssu];
W = [Wy;Wu];

% Objective J = [xvec' uvec']*daug(Qh,Rh)*[xvec;uvec]
%Bhi = inv(Ia)*Bh, A0i = inv(Ia)*A0,
%J = uvec'*Bhi'*Qh*Bhi*uvec + uvec'*Rh*uvec+2*x0'*A0i'*Qh*Bhi*uvec + (term in
%x0 only)

Bhi = Iai*Bh;
A0i = Iai*A0;

Hess = Bhi'*Qh*Bhi+Rh;
f0 = Bhi'*Qh*A0i; %To be post-multiplied by x0 at runtime

%Completing the square: 1/2*uvec'*Hess*uvec+f'*uvec =
%1/2*(uvec-hvec)'*Hess+(uvec-hvec) + (constant term); 
% hvec = -inv(Hess)*f = -inv(Hess)*f0*x0
% zvec = uvec-hvec;
%Gu \leq Sx0+W => Gz \leq Sx0+W-G*hvec

Wz = W;

Gz = G;

Sz = [Ss+G*inv(Hess)*f0];


% Assembling for use of upper and lower bounds on u by MOSEK quadprog.
high =hu(1:n*nu);
low =-hu(n*nu+1:2*n*nu);
% For mskqpopt (not very successful)
nyc = size(hy,1);
blc = ones(nyc,1).*(-inf);
