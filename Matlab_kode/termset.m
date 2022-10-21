function [Hn,hn] = termset(A,B,C,K,ymax,ymin,umax,umin)
%returnerer det settet vi har lyst til å ende opp inni, blir bare satt opp
%i MPC som ekstra beskrankninger som gjelder når N går mot uendelig
%the inputs are given as the discrete system dynamics:
%x(k+1)=Ax(k)+Bu(k)
%y(k)=Cx(k)
%K
%ymax, ymin, umax, umin
%Calculating terminal set for MPC
%Some of this could be done more elegantly with the 'polytope' routines in
%the MPT toolbox
%


nx = size(A,1); %Query only the length of the first dimension of A.
nu = size(B,2); %Query only the length of the second dimension of B.
ny = size(C,1); %Query only the length of the first dimension of C.

big = 1000;

H0 = [eye(nx);-eye(nx)];   %Defining constraints as H*x \leq h ??, making identity matrix 
h0 = big*ones(2*nx,1); %creating matrix with dimension 2*nx rows and 1 column containing 1000

Hi = [C;-C;-K;K];
hi = [ymax;-ymin;umax;-umin]; %constraints 
nci = size(hi,1);

ni = 10;
nstep = -1;
while (ni>0)
    ni = 0;
    nstep = nstep + 1;
    if (nstep == 0)
        Dyn = eye(nx); %first time the while loop runs, initialize an I-matrix
    else
        Dyn = (A-B*K)*Dyn; 
    end 
    for ik = 1:nci %iterating trough the constraints
        M= (Hi(ik,:)*Dyn)
        f = (Hi(ik,:)*Dyn)'; %creating objective function to minimize
        %How is it created?
        %[xopt,fval,lambda,exitflag,how]=mpt_solveLP(-f,H0,h0,[],[],[],0);
        %f=objective function as vector, 
        [x,fval,exitflag,output,lambda]=linprog(-f,H0,h0); %solves min -f'*x such that H0*x ≤ h0
        %x=optimal solution, fval = objective function value at minimum, lambda=lagrange multipliers
        if (-fval > hi(ik)) %one of the constraints is active
            ni = ni + 1; 
            H0 = [H0;f']; 
            l = hi(ik)
            h0 = [h0;hi(ik)];
        end
    end
end

%Removing redundant constraints from initial 'box'
nc0 = size(H0,1);

H1 = H0(2*nx+1:nc0,:);
h1 = h0(2*nx+1:nc0);

for ik = 1:2*nx
        f = H0(ik,:)';
        [xopt,fval,lambda,exitflag,how]=linprog(-f,H0,h0);%forskjell fra linprog?
        if (-fval > h0(ik))
            %nbox = nbox + 1
            H1 = [H1;f'];
            h1 = [h1;h0(ik)];
        end 
end
Hn = H1  %Could also have checked for redundancy of other constraints 
hn = h1




            