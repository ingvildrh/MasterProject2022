N = 1000;

G = -1 + (2).*rand(N,3);
Gn = normr(G); %normalizes a matrix 
Gt=Gn*Gn'; %
y = -10 +20*rand(N,1);
y0 = y;

Qinv = eye(N);
MatAcset = eye(N);
Aset = zeros(1,N);

%incompatible elements
ysign = sign(y)'; %returns same y but where all pos elements are 1, all negative are -1, 0 for 0
signtest=-(-1).^(Aset); %hvorfor i alle dager skal han ha den her?
posneg = ysign-signtest;%positive values correspond to candidates to enter the basis
indEnter2 = find(ysign>0); %returns all the indices of the positive elements in y 
indEnter = find(posneg>0);
indLeave = find(posneg<0);%returns all the indices of the negative elements in y
%Question: this won´t give the negative elements, since they are now 0 in
%posneg

ynorm = 0;

yabstt = sum(abs(y(union(indEnter,indLeave))));

%yd = norm(y-yopt);

yabs = [yabstt];
%this loop runs as long as we have elements that needs to leave or be added
%to the active set
while ~isempty(union(indEnter,indLeave)) %union(a,b) returns the data in both a and b with no repetitions
    disp('her')
    if ~isempty(indLeave) %if we have elements that needs to leave
        
        [v,j] = min(y(indLeave));
        i = indLeave(j);
        
        Mat = (MatAcset-Gt);
        Mvec = Mat(:,i);
        
        yplus = y - y(i)*(1+Mvec(i))^(-1)*Qinv*Mvec;
        Qinvplus = Qinv - (1+Mvec(i))^(-1)*Mvec*Qinv(i,:);
        Aset(i) = 0;
        MatAcset(i,i) = 1;
        
        
    elseif ~isempty(indEnter) %if we have elements that needs to be added
        %makes an element enter the base
        [v,j] = max(y(indEnter)); %finding the highest element v, and its index j
        i = indEnter(j);
        
        Mat = (MatAcset-Gt);
        Mvec = Mat(:,i);
        
        yplus = y - y(i)*(-1+Mvec(i))^(-1)*Qinv*Mvec; %equation 13, update y
        Qinvplus = Qinv - (-1+Mvec(i))^(-1)*Mvec*Qinv(i,:); %equation 12
        Aset(i) = 1;
        
        MatAcset(i,i) = 0;
        
    end



    %recomputes the incompatible elements
    ysign = sign(yplus)';%returns same y but where all pos elements are 1, all negative are -1, 0 for 0
    signtest=-(-1).^(Aset);
    posneg = ysign-signtest;%positive values correspond to candidates to enter the basis
    indEnter = find(posneg>0);
    indLeave = find(posneg<0);

    %test = union(indEnter,indLeave);

    yabsloop = sum(abs(y(union(indEnter,indLeave)))) %sum of all the indices in y that 
    ynormloop = norm(yplus-y)



    y = yplus;
    Qinv = Qinvplus;


    yabs = [yabs yabsloop];
    ynorm = [ynorm ynormloop];

    %%%distance to optimal
    %ydloop = norm(y-yopt)
    %yd = [yd ydloop];

end

plot(yabs)
%plot(ynorm)
%plot(yd)




