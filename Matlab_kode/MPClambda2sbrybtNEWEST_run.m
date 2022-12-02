tend = 100;
%  Using the matric inversion lemma for updating inv(Q)
%  Trying sparse implementation.
% Version b: Always delete inactve constraint before adding new active
% constraint.  Select largest violation of complementarity for
% removal/addition of constraints
% Version br: Removed unreliable LP formulation to check for infeasibility.
%  Rely solely on the size of qdiv instead.
% Version y: Update y as in (17) of draft paper
% Version bryb: minor tweaks
% Version brybt: Warmstarting from solution at previous timestep.
% Calculating times for optimization and for cleaning Q matrix separately
if (Model == 1)
    xinit = [5;-2];
    tosave1 = zeros(1,tend);
    while_iterations1 = zeros(1, 100);
elseif (Model == 2)
    xinit = 1.0*[25.5724;25.3546;9.7892;0.2448]; %(1.0: Feasible. 1.04 Feasible. 1.1: Infeasible, (but cycles indefinitely _No, not with updated limits on qdiv).  qdiv not very small. 
                                                 %1.5: Infeasibility caught through qdiv \approx 0. 
    while_iterations2 = zeros(1, 100);                                         
elseif (Model ==3)
  %  xinit = 0.1*randn(nx,1);
  load xinit3
  xinit = 4.98*xinit;   %Feasible.  One measurement violates constraints at t=0, but measurement constraints are only included from t=1.
  %xinit = 5.0*xinit; %Infeasible
  tosave2 = zeros(1,tend);
end

x0 = xinit;

xsave = zeros(nx,tend+1);
usave = zeros(nu,tend);
zsave = zeros(nu,tend);
lamsave = zeros(316, tend);
%tosave = zeros(1,tend);
%tcsave = zeros(1,tend);
feasflag = 1;  %feasflag = 0 indicates infeasible problem.

xsave(:,1)=xinit;

nc = size(Gz,1);
Hinv = inv(Hess);
%IGi = inv(eye(nc)-Gz*Hinv*Gz');
IGi = (eye(nc)-Gz*Hinv*Gz');
IGis = sparse(IGi);
hif = Hinv*f0;

 %  uvec = -Hinv*Gz'*lam - hif*x0;
  HGz = -Hinv*Gz';
  HGzu = HGz(1:nu,:);
  hifu = -hif(1:nu,:);

actset = sparse(logical(zeros(nc,1)));
actset0 = actset;
Qmat0i = speye(nc);

tt = tic;
%for ik = 1:1

tot = tic;
for ik = 1:tend
    %ik;
    %Implicit MPC problem.  Hoping to do this efficiently...
    
    %Starting with all constraints inactive
    to = tic;
    
    feasflag = false;
    solved = false;
    ix = 0;
    
    Qmat0i = speye(nc);
    actset = actset0;
   
    y0 = -Sz*x0-Wz;
    

    while (~solved)
        ix = ix +1;
        
        if (ix == 1)
            y = y0;
        else
            y = y0-y0(iz)*vAd;
            y0 = y;
        end
        
        lam = y.*actset;
        [i1,i1z] = min(lam);        
        if(i1 >= 0)
            i1 = [];
        end
        
        if(~isempty(i1))
            iz = i1z;
            actset(iz) = 0;
            qc = 1; % Remove constraint from active set
        else
            [i2,i2z] = max(y-lam);
            if (i2 <= 0)
                i2 = [];
            end
            if(~isempty(i2))
                iz = i2z;
                actset(iz) = 1;
                qc = -1; % Add constraint to active set
            else
                iz = [];
            end
        end       
        
        if(~isempty(iz))
                qu = IGis(:,iz);
                vA = Qmat0i*qu;
                qdiv = qc+vA(iz); %qc = 1/qc
 
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

                vAd = (1/qdiv)*vA;
                Qmat1i = Qmat0i-vAd*Qmat0i(iz,:);
                Qmat0i = Qmat1i;                   
        else
            solved = true;
            feasflag = true;

        end
        
    end
    tk = toc(to);
    if Model == 1
        while_iterations1(1,ik) = ix;
    end
    if Model == 2
        while_iterations2(1,ik) = ix;
    end
    if (feasflag == 0)
        break
    end
    
 %   lam = actset.*y; %Pick out the positive elements in y
    
  %   uvec = -Hinv*Gz'*lam - hif*x0;  
  %  u = uvec(1:nu);
    lamsave(:, ik+1)=lam;
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
    
  % tc = tic;
  % Qmat0i=clean(Qmat0i,1e-12);
   % tc = toc(tc);
   % tcsave(1,ik) = tc;
end
total_time = toc(tot);

toc(tt)


if (Model == 1)
    figure(1)
    tx = linspace(0,tend,tend+1);
    plot(tx,xsave(1,:),'b',tx,xsave(2,:),'r')
    legend('x1', 'x2')
    ylabel('x')
    xlabel('time step')
    title('State trajectories')

    figure(2)

    tu = linspace(0,tend-1,tend);
    plot(tu,usave)
    ylabel('u')
    xlabel('time step')
    title('System input')
    
    
elseif (Model == 2)
    ym = C*xsave;                        
    figure(1)
    tx = linspace(0,tend,tend+1);
    plot(tx,ym(1,:),'b',tx,ym(2,:),'r')
    legend('ym1', 'ym2')
    xlabel('time step')
    ylabel('ym')
    title('Output trajectories')
    grid


    figure(2)
    tu = linspace(0,tend-1,tend);
    plot(tu,usave(1,:),'b',tu,usave(2,:),'r')
    legend('u1', 'u2')
    ylabel('u')
    xlabel('time step')
    title('Input trajectories')
    grid
elseif (Model == 3)
    ym = C*xsave;                        
    figure(1)
    tx = linspace(0,tend,tend+1);
    plot(tx,ym(1,:),'b',tx,ym(2,:),'r',tx,ym(3,:),'k',tx,ym(4,:),'g')
    title('Output trajectories')

    figure(2)
    tu = linspace(0,tend-1,tend);
    plot(tu,usave(1,:),'b',tu,usave(2,:),'r',tu,usave(3,:),'k',tu,usave(4,:),'g')
end

    