%Version b.  Pulling out upper and lower bounds on u and giving these
%explicitly to qpOASES_sequence.

tend = 100;

if (Model == 1)
    xinit = [5;-2];
 %   xinit = [5;-4];
     load Model1feasSet
     xinit = 0.9999*P10V(10,:)';
elseif (Model == 2)
    xinit = 1.0*[25.5724;25.3546;9.7892;0.2448];
elseif (Model ==3)
  %  xinit = 0.1*randn(nx,1);
  load xinit3
  xinit = 4.98*xinit;
end

x0 = xinit;

xsave = zeros(nx,tend+1);
usave = zeros(nu,tend);
zsave = zeros(nu,tend);

xsave(:,1)=xinit;

f = zeros(n*nu,1);
hif = inv(Hess)*f0;

opt = qpOASES_options('MPC');
ncGy = size(Wy,1);
lb = -1e6*ones(ncGy,1);

tic

for ik = 1:tend
    %ik
   % if(ik==20)
   %     tic
   % end
   if (ik ==1)
        [QP,zvec,fval,exitflag,iter,lambda,auxOutput] = ...
                qpOASES_sequence( 'i',Hess,f0*x0,Gy,low,high,lb,Ssy*x0+Wy,opt);
 %   [zvec,fval,exitflag,output,lambda]=quadprog(Hess,f0*x0,Gy,Ssy*x0+Wy,[],[],low,high);
   else
        [zvec,fval,exitflag,iter,lambda,auxOutput] = ...
                       qpOASES_sequence( 'h',QP,f0*x0,low,high,lb,Ssy*x0+Wy,opt);
   end
    if (exitflag < 0)
        toc
        exitflag
        break
    end
   % if(ik==20)
   %     toc
   % end
    
    z = zvec(1:nu);
  %  u = z(1:nu)-hif(1:nu,:)*x0;
    u = z(1:nu);

    x1 = A*x0+B*u;
    
    x0 = x1;
    xsave(:,ik+1)=x0;
    usave(:,ik) = u;
    zsave(:,ik) = z;
end

toc

 qpOASES_sequence( 'c',QP )
 

if (Model == 1)
    figure(1)
    tx = linspace(0,tend,tend+1);
    plot(tx,xsave(1,:),'b',tx,xsave(2,:),'r')
    title('Output trajectories')

    figure(2)
    tu = linspace(0,tend-1,tend);
    plot(tu,usave)
    
elseif (Model == 2)
   % figure(1)
   %tx = linspace(0,tend,tend+1);
   % plot(tx,xsave(1,:),'b',tx,xsave(2,:),'r',tx,xsave(3,:),'k',tx,xsave(4,:),'g')
   %title('State trajectories')
    
    y = C*xsave;                        
    figure(1)
    tx = linspace(0,tend,tend+1);
    plot(tx,y(1,:),'b',tx,y(2,:),'r')
    title('Output trajectories')


    figure(2)
    tu = linspace(0,tend-1,tend);
    plot(tu,usave(1,:),'b',tu,usave(2,:),'r')
    
elseif (Model == 3)
    y = C*xsave;                        
    figure(1)
    tx = linspace(0,tend,tend+1);
    plot(tx,y(1,:),'b',tx,y(2,:),'r',tx,y(3,:),'k',tx,y(4,:),'g')
    title('Output trajectories')

    figure(2)
    tu = linspace(0,tend-1,tend);
    plot(tu,usave(1,:),'b',tu,usave(2,:),'r',tu,usave(3,:),'k',tu,usave(4,:),'g')
end
    
 
    



    