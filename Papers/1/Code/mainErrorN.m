
close all
clear all

order = 2;  % order of the spline basis

ntrial = 10;
np = zeros(length(ntrial),1);% vector of number of agents. For every entry np(i) an initial datum
         %from a certain distribution is drawn and the cloud of point is
         %generated
bdim = zeros(length(ntrial),1);% dimension of the function space
for i = 1:ntrial
    np(i) = i+2;
    bdim(i) = 3*np(i)-5;
end

L = zeros(length(np),1);
ndim = 2;   % dimensionality of agents
nfun = 8;   % choice of influence function (check file influence.m)
Mconstr = 100;

ErrorN = zeros(length(np),1);
Error2 = zeros(length(np),1);
Error2true = zeros(length(np),1);
ErrorInf = zeros(length(np),1);
ErrorTraj = zeros(length(np),1);

ti = 0;     %initial time
tf = 0.5;   %final time
dt = 0.001; %maximum time step

movebox = -5;%[-5,-1,3];
sizebox = 5;%[1,1,1];

for trial = 1:length(np)
    
    disp(trial)
    tic
    
    %% generate initial data and simulate
    
    flag = 1;
    while flag
        flag = 0;
        init = zeros(np(trial),ndim);
        for j = 1:np(trial)
            init(j,:) = sizebox*rand(ndim,1)+movebox;
        end
        init = reshape(init,[np(trial)*ndim 1]);
        h = @(t,x)cscdyn(t,x,np(trial),ndim,nfun);
        %[t,x] = RK4(init,ti,tf,dt,h); %solver for the ODE
        options = odeset('MaxStep',dt);
        [t,x] = ode113(h,[ti,tf],init,options); %solver for the ODE
        for j = 1:np(trial)
            if (norm(x(end,j:j+1)) > 1e+10)||(isnan(norm(x(end,j:j+1))))
                flag = 1;
            end
        end
    end
    
    toc
    
%     figure
%     hold on
%     for j = 1:np(trial)
%         plot(x(1:end,j),x(1:end,j+1));
%     end
%     hold off

        %% we modify the data in order to retrieve distances and approximated velocities
    
    z = x;
    nsteps = length(t);
    x = zeros(np(trial),ndim,nsteps);
    for s = 1:nsteps
        temp = squeeze(z(s,:));
        x(:,:,s) = reshape(temp,[np(trial) ndim]);
    end

    xdot = zeros(np(trial),ndim,nsteps); %velocities are approximated by finite differences
    for i = 1:np(trial)
        for m = 1:nsteps-1
            xdot(i,:,m) = (x(i,:,m+1)-x(i,:,m))/(t(m+1)-t(m));
        end
        xdot(i,:,nsteps) = xdot(i,:,nsteps-1);
    end


    X = zeros(np(trial),np(trial),ndim,nsteps);   %matrix of distance differences
    D = zeros(np(trial),np(trial),nsteps);        %matrix of distances between points
    for i = 1:np(trial)
        for j = 1:np(trial)
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
    %x = z;

    A = squeeze(max(D(:,:,1:nsteps)));
    L(trial) = max(squeeze(max(A(:,1:nsteps))));   %maximum x-value displayed

    %array with distances of agents sorted in increasing order
    %sortedD = sort(reshape(D,[np(trial)*np(trial)*nsteps 1]));
    %sortedD = sortedD(sortedD > 0);
    
    %% THE PROCEDURE
    
    %knot vector for the spline basis: it is evenly distributed. The
    %first and last n = order knots are 0, to have a jump of dicontinuity
    %at 0 and at L. we add an external node since it is unknown why the
    %last approximation function has always coeff = 0

    stepKnot = L(trial)/(bdim(trial) + order - 4);
    basicKnots = 0:stepKnot:L(trial);
    knots_0 = [0,0,basicKnots,L(trial),L(trial)];

    toc

    [d,C] = fitDynamics(xdot,X,D,bdim(trial),knots_0);
    [b,h,p] = size(C);
    d = reshape(d,1,h*p);
    C = reshape(C,bdim(trial),h*p);
    
    diff = zeros(bdim(trial),bdim(trial));
    for con = 1:bdim(trial)
        if con < bdim(trial)
            diff(con,con) = 1;
            diff(con,con+1) = -1;
        end
    end

%     least square approximation
%     y_0(:,trial) = lsqlin(C',d',[],[]);
    cvx_begin
        cvx_precision best
        variable x(bdim(trial),1)
        minimize( norm(C'*x - d',2) )
        subject to
            2*norm(x,Inf) + norm(diff*x,Inf) <= Mconstr;
    cvx_end
    y_0 = x;
    
    ErrorN(trial) = 1/length(t)*1/np(trial)*norm(C'*x - d',2)^2;
    
    close all
    dd = 0:0.001:L(trial);
    bb = influence(dd,nfun);
    aa = aLearned(bdim(trial),knots_0,y_0,dd);
    %hold on
    %plot(dd,bb,'r');
    %plot(dd,aa,'b');
    %hold off
    %drawnow
%     for index = 1:length(dd)
%         Error2(trial) = Error2(trial) + norm(aa(i) - bb(i),2)*dd(i)^2;
%     end
    
    D = reshape(D,[np(trial)*np(trial)*length(t) 1]);
    bb = influence(D,nfun);
    aa = aLearned(bdim(trial),knots_0,y_0,D);
    for index = 1:length(D)
        Error2true(trial) = Error2true(trial) + (norm(aa(i) - bb(i),2)*D(i))^2;
        %Error2(trial) = Error2(trial) + (norm(aa(i) - bb(i),2))^2;
    end
    Error2true(trial) = 1/length(t)*1/np(trial)*1/np(trial)*Error2true(trial)
    %Error2(trial) = 1/tf*1/np(trial)*1/np(trial)*Error2(trial)
    
    %ErrorInf(trial) = max(abs(aa - bb));
    
    toc
    
%     h = @(t,x)cscdynApprox(t,x,np(trial),ndim,bdim(trial),knots_0,y_0);
%     [t,x] = RK4(init,ti,tf,dt,h); %solver for the ODE
%     
%     cc = zeros(length(t),1);
%     for s = 1:length(t)
%         for j = 1:np(trial)
%             cc(s) = cc(s) + 1/np(trial)*norm(x(s,j:j+1) - z(s,j:j+1),1);
%         end
%     end
%     ErrorTraj(trial) = max(cc);
%     
%     toc
    
%     figure
%     hold on
%     for j = 1:np(trial)
%         plot(x(1:end,j),x(1:end,j+1));
%     end
%     hold off

%     
%     figure(trial)
%     hold on
%     plot(dd,bb,'r');
%     plot(dd,aa,'b');
%     hold off
    
    toc
    
end


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
