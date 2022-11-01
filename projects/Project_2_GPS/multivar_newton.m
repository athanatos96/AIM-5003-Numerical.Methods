% Multivariate Newton's Method for GPS project.
% 
% Description
% -----------
% Solve the GPS nonlinear system of equations 
% for the unknown (x,y,z,d) using Newton's method.
%
% \param[in]  x0 Initial vector
% \param[in]  S1 Position of first satellite
% \param[in]  S2 Position of second satellite
% \param[in]  S3 Position of third satellite
% \param[in]  S4 Position of fourth satellite
% \param[in]  T  Measured time intervals
% \param[in]  k  Maximum number of iterations of the method.
%
% \param[out] x  Approximate solution.

function x = multivar_newton(x0,S1,S2,S3,S4,T,k)
x=x0;
% loop for Newton's iteration
for j=1:k
  F=GPS_f(x,S1,S2,S3,S4,T);   % Function at x
  DF=GPS_df(x,S1,S2,S3,S4,T); % Jacobian at x
  
  % Check that Jacobian matrix is invertible.
  if(rank(DF)<length(x0))
    disp("Jacobian is singular! Exiting.\n")
    return
  endif
  
  % Instead of computing the inverse of the Jacobian, solve the equivalent 
  % linear system
  s=DF\F;
  x=x-s;
  printf("x=(%f, %f, %f, %f)\n", x(1), x(2), x(3), x(4));
end % End of iteration loop