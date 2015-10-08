tic
np = 50; % vector of number of agents. For every entry np(i) an initial datum
        %from a certain distribution is drawn and the cloud of point is
        %generated
ndim = 2;   % dimensionality of agents
nfun = 8;  % choice of influence function (check file influence.m)

ti = 0;   %initial time
tf = 0.5;  %final time
dt = 0.01; %maximum time step

movebox = -5;%[-5,-1,3];
sizebox = 15;%[2,2,2];

x00=cell(length(np));
for i = 1:length(np)
    init = zeros(np(i),ndim);
    for j = 1:np(i)
        init(j,:) = sizebox*rand(ndim,1)+movebox;
%         R = unidrnd(3,1,1);
%         if  R == 1
%             init(j,:) = sizebox(1)*rand(ndim,1)+movebox(1); %each entry uniformly
%                                                 %taken from [movebox,sizebox+movebox]
%         elseif R == 2
%             init(j,:) = sizebox(2)*rand(ndim,1)+movebox(2);
%         else
%             init(j,:) = sizebox(3)*rand(ndim,1)+movebox(3);
%         end
    end
    x00{i} = init;
end

for i = 1:length(np)
    init = x00{i};
    init = reshape(init,[np(i)*ndim 1]);
    h = @(t,x)cscdyn(t,x,np(i),ndim,nfun);
    options = odeset('MaxStep',dt);
    [t,x] = ode113(h,[ti,tf],init,options); %solver for the ODE

    %save the relevant data into cloud.mat
    if length(np) == 1
        fullFileName = sprintf('cloud.mat');
    else
        fullFileName = sprintf('cloud%d.mat',i);
    end
    tempnp = np;
    np = np(i);
    save(fullFileName,'t','x','np','ndim','nfun','ti','tf','dt');
    np = tempnp;
end

toc
