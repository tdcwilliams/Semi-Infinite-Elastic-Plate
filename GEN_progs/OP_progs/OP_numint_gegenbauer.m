function [t,w]=OP_numint_gegenbauer(alpC,Nint);

[t,w]=OP_numint_jacobi(alpC-.5,alpC-.5,Nint);