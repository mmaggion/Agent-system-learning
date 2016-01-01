function [x]=cscdynApprox(time,xinit,np,ndim,dBasis,knots,coeff)

a=zeros(np,np,ndim);
distance=zeros(np,np);
difpos=zeros(np,np,ndim);
x0=reshape(xinit,[np ndim]);
x0=x0(1:np,:);

for i=1:np
    difpos(i,:,:)=x0-repmat(x0(i,:),np,1);
    distance(i,:)=sqrt(sum((repmat(x0(i,:),np,1)-x0).^2,2));
end

distance(logical(eye(size(distance)))) = 0.01; %avoids division by zero

for k=1:ndim
    weights = aLearned(dBasis,knots,coeff,reshape(distance,[np*np 1]));
    %a(:,:,k) = 1/np*reshape(weights,[np np]).*difpos(:,:,k)./distance(:,:);
    a(:,:,k) = 1/np*reshape(weights,[np np]).*difpos(:,:,k);
end

x(1:np,:)=squeeze(sum(a,2));
x = reshape(x,[np*ndim 1]);
end 