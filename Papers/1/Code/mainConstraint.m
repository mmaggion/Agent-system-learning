
close all

order = 2;  % order of the spline basis

np = 50; % vector of number of agents. For every entry np(i) an initial datum
         %from a certain distribution is drawn and the cloud of point is
         %generated
bdim = 2*np;   % dimension of the function space
         
%y_0 = zeros(bdim,length(np));
%knots_0 = zeros(length(np),bdim + order);
ndim = 2;   % dimensionality of agents
nfun = 8;   % choice of influence function (check file influence.m)
Mconstr = 0:10:200;
ErrorM = zeros(length(Mconstr),1);

ti = 0;     %initial time
tf = 0.5;   %final time
dt = 0.001; %maximum time step

movebox = -3;%[-5,-1,3];
sizebox = 6;%[1,1,1];

flag = 1;
while flag
    flag = 0;
    init = zeros(np,ndim);
    for j = 1:np
        init(j,:) = sizebox*rand(ndim,1)+movebox;
    end
    init = reshape(init,[np*ndim 1]);
    h = @(t,x)cscdyn(t,x,np,ndim,nfun);
    options = odeset('MaxStep',dt);
    [t,x] = ode113(h,[ti,tf],init,options); %solver for the ODE
    for j = 1:np
        if norm(x(end,j:j+1)) > 1e+10
            flag = 1;
        end
    end
end

    %% we modify the data in order to retrieve distances and approximated velocities

z = x;
nsteps = length(t);
x = zeros(np,ndim,nsteps);
for s = 1:nsteps
    temp = squeeze(z(s,:));
    x(:,:,s) = reshape(temp,[np ndim]);
end

xdot = zeros(np,ndim,nsteps); %velocities are approximated by finite differences
for i = 1:np
    for m = 1:nsteps-1
        xdot(i,:,m) = (x(i,:,m+1)-x(i,:,m))/(t(m+1)-t(m));
    end
    xdot(i,:,nsteps) = xdot(i,:,nsteps-1);
end


X = zeros(np,np,ndim,nsteps);   %matrix of distance differences
D = zeros(np,np,nsteps);        %matrix of distances between points
for i = 1:np
    for j = 1:np
        for m = 1:nsteps
            D(i,j,m) = norm(x(i,:,m) - x(j,:,m),2);
            if i == j
                X(i,j,:,m) = 0;
            else
                %X(i,j,:,m) = (x(i,:,m) - x(j,:,m))/D(i,j,m);
                X(i,j,:,m) = x(i,:,m) - x(j,:,m);
            end
        end
    end
end
x = z;

A = squeeze(max(D(:,:,1:nsteps)));
L = max(squeeze(max(A(:,1:nsteps))));   %maximum x-value displayed

stepKnot = L/(bdim - order - 1);
basicKnots = 0:stepKnot:L;
knots_0 = [0,0,basicKnots,L,L];


[d,C] = fitDynamics(xdot,X,D,bdim,knots_0);
[b,h,p] = size(C);
d = reshape(d,1,h*p);
C = reshape(C,bdim,h*p);

diff = zeros(bdim,bdim);
for con = 1:bdim
    if con < bdim
        diff(con,con) = 1;
        diff(con,con+1) = -1;
    end
end

for trial = 1:length(Mconstr)
    
    disp(trial)
    tic
    
    %% generate initial data and simulate
    
    

    %array with distances of agents sorted in increasing order
    %sortedD = sort(reshape(D,[np(trial)*np(trial)*nsteps 1]));
    %sortedD = sortedD(sortedD > 0);
    
    %% THE PROCEDURE
    
    %knot vector for the spline basis: it is evenly distributed. The
    %first and last n = order knots are 0, to have a jump of dicontinuity
    %at 0 and at L. we add an external node since it is unknown why the
    %last approximation function has always coeff = 0

    %least square approximation
    %y_0(:,trial) = lsqlin(C',d',[],[]);
    cvx_begin
        cvx_precision best
        variable x(bdim,1)
        minimize( norm(C'*x - d',2) )
        subject to
            2*norm(x,Inf) + norm(diff*x,Inf) <= Mconstr(trial);
    cvx_end
    
    ErrorM(trial) = 1/length(t)*1/np*norm(C'*x - d',2)
    
    dd = 0:0.001:L;
    bb = influence(dd,nfun);
    aa = aLearned(bdim,knots_0,x,dd);
    
    figure(trial)
    hold on
    plot(dd,bb,'r');
    plot(dd,aa,'b');
    hold off
    drawnow
    
    toc
    
end

plot(Mconstr,ErrorM);


% dd = 0:0.001:min(L);                                %approximation interval
% bb = influence(dd,nfun);                            %original function
% aa = zeros(length(np),length(dd));
% figure
% hold on
% for trial = 1:length(np)
%     aa(trial,:) = aLearned(bdim,knots_0(trial,:),y_0(:,trial),dd);    %reconstructed functions
% end
% muhat = zeros(1,length(dd));
% sigmahat = zeros(1,length(dd));
% muci = zeros(2,length(dd));
% sigmaci = zeros(2,length(dd));
% for step = 1:length(dd)
%     [muhat(step),sigmahat(step),muci(:,step),sigmaci(:,step)] = normfit(aa(:,step));
% end
% plot(dd,bb,'g',...
%      dd,muhat,'b',...
%      dd,muci(1,:),'r:',...
%      dd,muci(2,:),'r:');
% hold off
