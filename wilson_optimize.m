%wilson_optimize   Find most ill conditioned Wilson-like matrices.

n_min = 1;
n_max = 10;
tol = 1e-14; % For testing singularity and definiteness.

cond_max = 0; ind_max = 0;
cond_max_pd = 0; ind_max_pd = 0;
n_sing = 0; n_num_sing = 0;

A = zeros(4);
for a = n_min:n_max
  A(1,1) = a;
  for b = n_min:n_max
    A(1,2) = b;
    A(2,1) = b;
    for c = n_min:n_max
      A(1,3) = c;
      A(3,1) = c;
      for d = n_min:n_max
        A(1,4) = d;
        A(4,1) = d;
        for e = a:n_max
          A(2,2) = e;
          for f = n_min:n_max
            A(2,3) = f;
            A(3,2) = f;
            for g = n_min:n_max
              A(2,4) = g;
              A(4,2) = g;
              for h = e:n_max
                A(3,3) = h;
                for i = n_min:n_max
                  A(3,4) = i;
                  A(4,3) = i;
                  for j = h:n_max
                    A(4,4) = j;
                    detA = d^2*f^2 - 2*c*d*f*g + c^2*g^2 -...
                           d^2*e*h + 2*b*d*g*h - a*g^2*h +...
                           2*c*d*e*i - 2*b*d*f*i -...
                           2*b*c*g*i + 2*a*f*g*i + b^2*i^2 -...
                           a*e*i^2 - c^2*e*j + 2*b*c*f*j -...
                           a*f^2*j - b^2*h*j + a*e*h*j;

                    if detA == 0
                        % fprintf('%11.0f singular\n',i)
                       n_sing = n_sing + 1;
                       continue
                    end
                    evals = eig(A);
                    s = abs(evals);  % Singular values.
                    if min(s) < tol
                        % fprintf('%11.0f "numerically" singular\n',i)
                        n_num_sing = n_num_sing + 1;
                        continue
                    end
                    kappa = max(s)/min(s);
                    if kappa > cond_max
                       A_max = A;
                       cond_max = kappa;
                       fprintf(['a = %2.0f, b = %2.0f, c = %2.0f, d = %2.0f' ...
                                '  cond_max = %9.2e\n'],a,b,c,d,cond_max)
                    end
                    if kappa > cond_max_pd
                       if min(evals) >= tol
                         A_max_pd = A;
                         cond_max_pd = kappa;
                         fprintf(['a = %2.0f, b = %2.0f, c = %2.0f, d = %2.0f' ...
                                  '  cond_max_pd = %9.2e\n'],a,b,c,d,cond_max_pd)
                       end
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end

n_sing
n_num_sing
A_max
cond_max
A_max_pd
cond_max_pd
