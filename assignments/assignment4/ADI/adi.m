   clear; close all;

   a = 0; b=1;  c=0; d=1; n = 32;  tfinal = 0.5;

   m = n;
   h = (b-a)/n;        dt=h;   
   h1 = h*h;
   x=a:h:b; y=c:h:d;
   [X,Y] = meshgrid(x,y);

%-- Initial condition:

   t = 0;
   for i=1:m+1,
      for j=1:m+1,
         u1(i,j) = uexact(t,x(i),y(j));
      end
   end

%---------- Big loop for time t --------------------------------------

k_t = fix(tfinal/dt);

for k=1:k_t

t1 = t + dt; t2 = t + dt/2;

%--- sweep in x-direction --------------------------------------
tic;

for i=1:m+1,                              % Boundary condition.
  u2(i,1) = uexact(t2,x(i),y(1));
  u2(i,n+1) = uexact(t2,x(i),y(n+1));
  u2(1,i) = uexact(t2,x(1),y(i));
  u2(m+1,i) = uexact(t2,x(m+1),y(i));
end

for j = 2:n,                             % Look for fixed y(j) 

   A = sparse(m-1,m-1); b=zeros(m-1,1);
   for i=2:m,
      b(i-1) = (u1(i,j-1) -2*u1(i,j) + u1(i,j+1))/h1 + ...
		f(t2,x(i),y(j)) + 2*u1(i,j)/dt;
      if i == 2
        b(i-1) = b(i-1) + uexact(t2,x(i-1),y(j))/h1;
        A(i-1,i) = -1/h1;
      else
	if i==m
          b(i-1) = b(i-1) + uexact(t2,x(i+1),y(j))/h1;
          A(i-1,i-2) =  -1/h1;
	else
	   A(i-1,i) = -1/h1;
	   A(i-1,i-2) = -1/h1;
        end
      end

      A(i-1,i-1) = 2/dt + 2/h1;
    end

     ut = A\b;                          % Solve the diagonal matrix.
     for i=1:m-1,
	u2(i+1,j) = ut(i);
     end

 end                                    % Finish x-sweep.

%-------------- loop in y -direction --------------------------------

for i=1:m+1,                                % Boundary condition
  u1(i,1) = uexact(t1,x(i),y(1));
  u1(i,n+1) = uexact(t1,x(i),y(m+1));
  u1(1,i) = uexact(t1,x(1),y(i));
  u1(m+1,i) = uexact(t1,x(m+1),y(i));
end

for i = 2:m,

   A = sparse(m-1,m-1); b=zeros(m-1,1);
   for j=2:n,
      b(j-1) = (u2(i-1,j) -2*u2(i,j) + u2(i+1,j))/h1 + ...
                f(t2,x(i),y(j)) + 2*u2(i,j)/dt;
      if j == 2
        b(j-1) = b(j-1) + uexact(t1,x(i),y(j-1))/h1;
        A(j-1,j) = -1/h1;
      else
        if j==n
          b(j-1) = b(j-1) + uexact(t1,x(i),y(j+1))/h1;
          A(j-1,j-2) =  -1/h1;
        else
           A(j-1,j) = -1/h1;
           A(j-1,j-2) = -1/h1;
        end
      end

      A(j-1,j-1) = 2/dt + 2/h1;              % Solve the system
    end

     ut = A\b;
     for j=1:n-1,
        u1(i,j+1) = ut(j);
     end

 end                             % Finish y-sweep.

 t = t + dt;

%--- finish ADI method at this time level, go to the next time level.
      
end       %-- Finished with the loop in time
cpu_time = toc;
%----------- Data analysis ----------------------------------

  for i=1:m+1,
    for j=1:n+1,
       ue(i,j) = uexact(tfinal,x(i),y(j));
    end
  end

  e = max(max(abs(u1-ue)));        % The infinity error.
  fprintf('Max error of ADI with m=n=%d: %.15f\n',m,e);
  fprintf('CPU time of ADI with m=n=%d: %.6f seconds\n',m,cpu_time);

  figure(1)
  hold on;
  subplot(1,2,1);
  mesh(X,Y,u1');                    % Plot the computed solution.
  xlabel('x'); ylabel('y'); zlabel('u(x,y,t)');
  title(['Numerical solution by ADI Method'...
      ' with m=n=' num2str(m) ' grid intervals'])

  subplot(1,2,2);
  mesh(X,Y,abs(u1'-ue'));
  xlabel('x'); ylabel('y'); zlabel('Error');% Mesh plot of the error 
  title(['Error Distribution of ADI Method ' ...
       'with m=n=' num2str(m) ' grid intervals']);