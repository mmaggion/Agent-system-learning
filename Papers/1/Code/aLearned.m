function a = aLearned(basisD,knotVector,c,x)
%aLearned returns the value at x of the linear combination of spline basis functions
%basisD is the dimension of the function space
%knotVector is the vectors of knots used to create the spline basis
%c is the vector of coefficients for the linear combination
%x is the point at which we evaluate functions

a = 0;
orderB = length(knotVector)-basisD; %recover desidered order of splines

for i = 0:length(c)-1
    
    a = a + c(i+1)*bspline_basis(i,orderB,knotVector,x);

end