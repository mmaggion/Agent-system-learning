function [t,x]=RK4(x0,T0,Tend,Tstep,fun)
%% solving x'=f(t,x) with RK4
%% x0 line vector of size dim
%% x matrix (N+1)*dim


t = [T0:Tstep:T0+Tend];
x(1,:) = x0';
for k = 1:length(t)-1
   p1 = feval(fun,t(k),x(k,:)');
   p2 = feval(fun,t(k) + Tstep/2,x(k,:)' + Tstep/2*p1);
   p3 = feval(fun,t(k) + Tstep/2,x(k,:)' + Tstep/2*p2);
   p4 = feval(fun,t(k+1),x(k,:)' + Tstep*p3);
   x(k+1,:) = x(k,:) + Tstep*(1/6*p1 + 1/3*p2 + 1/3*p3 + 1/6*p4)';
end
