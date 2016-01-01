close all
clear all

order = 2;  % order of the spline basis


np = 20;% vector of number of agents. For every entry np(i) an initial datum
         %from a certain distribution is drawn and the cloud of point is
         %generated
bdim = 150;% dimension of the function space

ndim = 2;   % dimensionality of agents
nfun = 6;   % choice of influence function (check file influence.m)
Mconstr = 1/8*[10,15,20,25,30,35,40];
Xtraj = cell(length(Mconstr),1);

ti = 0;     %initial time
tf = 1;%0.5;   %final time
dt = 0.01; %maximum time step

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
    [t,x] = RK4(init,ti,tf,dt,h); %solver for the ODE
    %options = odeset('MaxStep',dt);
    %[t,x] = ode113(h,[ti,tf],init,options); %solver for the ODE
    for j = 1:np
        if (norm(x(end,j:j+1)) > 1e+10)||(isnan(norm(x(end,j:j+1))))
            flag = 1;
        end
    end
end

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

A = squeeze(max(D(:,:,1:nsteps)));
L = max(squeeze(max(A(:,1:nsteps))));
dd = 0:0.001:L;
bb = influence(dd,nfun);
aa = zeros(length(Mconstr),length(dd));

stepKnot = L/(bdim + order - 4);
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
    

    %% we modify the data in order to retrieve distances and approximated velocities
    
   %maximum x-value displayed

    %array with distances of agents sorted in increasing order
    %sortedD = sort(reshape(D,[np(trial)*np(trial)*nsteps 1]));
    %sortedD = sortedD(sortedD > 0);
    
    %% THE PROCEDURE
    
    %knot vector for the spline basis: it is evenly distributed. The
    %first and last n = order knots are 0, to have a jump of dicontinuity
    %at 0 and at L. we add an external node since it is unknown why the
    %last approximation function has always coeff = 0

    cvx_begin
        cvx_precision best
        variable x(bdim,1)
        minimize( norm(C'*x - d',2) )
        subject to
            2*norm(x,Inf) + norm(diff*x,Inf) <= Mconstr(trial);
    cvx_end
    y_0 = x;

    aa(trial,:) = aLearned(bdim,knots_0,y_0,dd);
    
    h = @(t,x)cscdynApprox(t,x,np,ndim,bdim,knots_0,y_0);
    [t,x] = RK4(init,ti,tf,dt,h); %solver for the ODE
    Xtraj{trial} = x;
    
    toc
    
end

C = autumn(length(Mconstr));

figure(1)
hold on
set(gca, 'color', [0.25 0.25 0.25])
for trial = 1:length(Mconstr);
    plot(dd,-aa(trial,:),'color',C(trial,:),'LineStyle','-.');
end
plot(dd,-bb,'w','LineWidth',1);
hold off

figure(2)
hold on
set(gca, 'color', [0.25 0.25 0.25])
for trial = 1:length(Mconstr);
    x = Xtraj{trial};
    for j = 1:np
        plot(x(1:end,j),x(1:end,j+1),'color',C(trial,:),'LineStyle','-.');
    end
end

for j = 1:np
    %plot(z(1:end,j),z(1:end,j+1),'color',C{1});
    plot(z(1:end,j),z(1:end,j+1),'w','LineWidth',1);
end
hold off