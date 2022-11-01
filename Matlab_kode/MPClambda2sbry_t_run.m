tend = 100; %number of steps to simulate 
%  Using the matric inversion lemma for updating inv(Q)
%  Trying sparse implementation.
% Version b: Always delete inactve constraint before adding new active
% constraint.  Select largest violation of complementarity for
% removal/addition of constraints
% Version br: Removed unreliable LP formulation to check for infeasibility.
%  Rely solely on the size of qdiv instead.
% Version y: Update y as in (17) of draft paper
%Version _t: adding timing info to version bry
if (Model == 1)
    xinit = [5;-2];
  %  load Model1feasSet
  %  xinit = 0.9999*P10V(10,:)';
    xinit = [-6;1.1];
    tosave1 = zeros(1,tend);
elseif (Model == 2)
    xinit = 1.0*[25.5724;25.3546;9.7892;0.2448]; %(1.0: Feasible. 1.04 Feasible. 1.1: Infeasible, (but cycles indefinitely _No, not with updated limits on qdiv).  qdiv not very small. 
                                                 %1.5: Infeasibility caught through qdiv \approx 0. 
    tosave2 = zeros(1,tend);
elseif (Model ==3)
  %  xinit = 0.1*randn(nx,1);
  load xinit3
  xinit = 4.98*xinit;   %Feasible.  One measurement violates constraints at t=0, but measurement constraints are only included from t=1.
  %xinit = 5.0*xinit; %Infeasible
 % load xinits
 % xinit = 0.5*xinits; %0.501 infeasible
end


x0 = xinit;


xsave = zeros(nx,tend+1);%allocating space for the variables in each step
usave = zeros(nu,tend);
zsave = zeros(nu,tend);


feasflag = 1;  %feasflag = 0 indicates infeasible problem. Changes if it is not possible to find a solution

xsave(:,1)=xinit;

f = zeros(n*nu,1);


nc= size(Gz,1);
Hinv = inv(Hess);
%IGi = inv(eye(nc)-Gz*Hinv*Gz');
IGi = (eye(nc)-Gz*Hinv*Gz'); %(I-GH^(-1)G^(T))
IGis = sparse(IGi); %saving as a glissen matrix
hif = Hinv*f0; %-H^(-1)*F(T)^
%7a og 7b

%  uvec = -Hinv*Gz'*lam - hif*x0;
HGz = -Hinv*Gz'; %
HGzu = HGz(1:nu,:); 
hifu = -hif(1:nu,:);

actsets = sparse(zeros(nc,1)); %allocating space for active set
Qsp0 = speye(nc); %speye returns a sparse nbyn identity matrix with ones on diagonal


%tt=tic;
%for ik = 1:1
%algorithm 1 solving for tend timesteps
while_iterations = zeros(1, 100);

for ik = 1:tend
    %ik;
    %Implicit MPC problem.  Hoping to do this efficiently...
    
    %Starting with all constraints inactive
    to=tic; %to measure elapsed time
    
    feasflag = false;

    actset = actsets;  %Initialized to all zeros. Why do we empty it every step?
    solved = false;
    izold = 0;
    ix = 0;
    
    Qmat0i = Qsp0;
    y0 = -Sz*x0-Wz; %equation 9?
    num_iterations = 0;
    while (~solved)
        num_iterations = num_iterations + 1;
        ix = ix +1;
        
        if (ix == 1) %first time, make sure to set y to y0
            y = y0;
        else %not first time, update y (find the eq)
      
            y = y0-y0(iz)*vAd;
            y0 = y;
        end %QP løsningen er i denne while loopen
        
        %[i1,i1z] = min(y.*actset);
        lam = y.*actset; %hvis det aktive settet er riktig blir dette lambda
        
        [i1,i1z] = min(lam); %returns iz as the lowest value in lam and i1z as the index of this value
        
        if(i1 >= -10*eps) %if the smallest value is above a certain limit, we remove it
            i1 = []; %is it sufficient to just think that i2<0?
        end
        
        if(~isempty(i1))
            iz  = i1z; %set the iz to the index
            actset(iz) = 0; %remove from the active set on that iz index
            qc = 1; % Remove constraint from active set
        else
             [i2,i2z] = max(y-lam); %is it sufficient to just think that i2>0?
           % [i2,i2z] = max(y.*(ones(nc,1)-actset));
            if (i2 <= 10*eps)
                i2 = [];
            end
            if(~isempty(i2)) %sett bp her
         %       iz  = i2;
                iz  = i2z;
                actset(iz) = 1; %the active set is binary
                qc = -1; % Add constraint to active set
            else  %(Only if i1=[] and i2 = []
                iz = [];
            end
        end

        
        if(~isempty(iz))
                qu = IGis(:,iz);
                vA = Qmat0i*qu;
                qdiv = qc+vA(iz); %qc = 1/qc
                vAd = (1/qdiv)*vA;
 
                %noe som ikke er nødvendig, ikke forstå
                if ((abs(qdiv) < 1e-13) | (abs(qdiv) > 1e12))
                    toc
                %Will only enter here (some times) when adding a new constraint    
                % Heuristically seems like a reliable test for infeasibility

                       feasflag = 0;
                       disp('Infeasible problem detected')

                       qdiv
                       feasflag = false;
                       break
                end
                S = vAd*Qmat0i(iz,:);
                o=vAd*Qmat0i(iz,:)
                Qmat1i = Qmat0i+vAd*Qmat0i(iz,:);
                Qmat0i = Qmat1i;                   
        else
            solved = true;
            feasflag = true;

        end
        
    end
    while_iterations(1,ik) = num_iterations;
    tk = toc(to);
    if (feasflag == 0)
        break
    end
    
    lam = actset.*y; %Pick out the positive elements in y
    
 %   uvec = -Hinv*Gz'*lam - hif*x0;  
  %  u = uvec(1:nu);
    
  %finne det første pådraget
    u = HGzu*lam+hifu*x0;
        
    % End solve implicit MPC problem

    x1 = A*x0+B*u;
    
    x0 = x1;
    xsave(:,ik+1)=x0;
    usave(:,ik) = u;
    if (Model ==1)
        tosave1(1,ik) = tk;
    end
    if (Model ==2)
        tosave2(1,ik) = tk;
    end
end
toc(tt)

if (Model == 1)
    figure(1)
    tx = linspace(0,tend,tend+1);
    plot(tx,xsave(1,:),'b',tx,xsave(2,:),'r')
    ylabel('x1')
    xlabel('timestep')
    title('State trajectories')

    figure(2)
    tu = linspace(0,tend-1,tend);
    plot(tu,usave)
    
elseif (Model == 2)
    y = C*xsave;                        
    figure(1)
    tx = linspace(0,tend,tend+1);
    plot(tx,y(1,:),'b',tx,y(2,:),'r')
    legend('y1', 'y2')
    title('Output trajectories')
    ylabel('y')
    xlabel('timestep')


    figure(2)
    ylabel('y')
    tu = linspace(0,tend-1,tend);
    plot(tu,usave(1,:),'b',tu,usave(2,:),'r')
    title('Input trajectories')
    legend('u1', 'u2')
    xlabel('timestep')
elseif (Model == 3)
    y = C*xsave;                        
    figure(1)
    tx = linspace(0,tend,tend+1);
    plot(tx,y(1,:),'b',tx,y(2,:),'r',tx,y(3,:),'k',tx,y(4,:),'g')
    legend('y1', 'y2')
    title('Output trajectories')
    xlabel('timestep')
    figure(2)
    tu = linspace(0,tend-1,tend);
    plot(tu,usave(1,:),'b',tu,usave(2,:),'r',tu,usave(3,:),'k',tu,usave(4,:),'g')
end

    