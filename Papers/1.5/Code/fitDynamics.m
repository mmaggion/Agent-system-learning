function [VDOT,DYN] = fitDynamics(derivatives,pointers,distances,basisD,knotVector)

%derivatives is a np*ndim matrix containing the accelereation of the system
%pointers is a np*np*ndim matrix containing the differences between the
%velocities of the particles
%distances is a np*npmatrix containing the distances between pairs of points
%basisD is the dimension of the function space
%knotVector is the vector of knots necessary for constructing a spline
%basis

[N,spaceD,T] = size(derivatives);
orderB = length(knotVector)-basisD;

VDOT = zeros(size(derivatives,1),size(derivatives,2)*size(derivatives,3));         %vector of enqueued derivatives
PDYN = cell(3,1);               %cell array to be then stuffed into DYN
DYN = zeros(basisD,N,spaceD*T); %vector of enqueued dynamics

idx = 1;
for m = 1:T
    for d = 1:spaceD
        VDOT(:,idx) = derivatives(:,d,m); %enqueue derivatives vectors
        idx = idx + 1;
    end
end

TDYN = zeros(N,T*spaceD);
for k = 0:basisD-1
    idx  = 1;
    for m = 1:T
        B = reshape(bspline_basis(k,orderB,knotVector,reshape(distances(:,:,m),[N*N 1])),[N N]);
        for d = 1:spaceD
            TDYN(:,idx) = 1/N*diag(B*pointers(:,:,d,m)); %enqueue dynamics vectors
            idx = idx + 1;
        end   
    end
    PDYN{k+1} = TDYN; %store dynamics vector in a cell array according to basis index
end

DYN = zeros(basisD,size(PDYN{1},1),size(PDYN{1},2));
for k = 1:basisD
    DYN(k,:,:) = PDYN{k}; %stack PDYN matrices according to basis index
end

end