clear; clc
N = 20; dt = 1/N; drind = (1.5 - (-2))/N;
pedelta = 4*sqrt(3)*sqrt(pi);

% Generate s values
for i = 0:N
    % As the number of segments increases, s should be appropriately adjusted
    % to reduce the singularity of the matrix. However, note that the overall
    % distribution range of s should remain close to [0, 1].
    s(i+1) = i/N;
    if s(i+1) > 1
        i  % Display index if s exceeds 1
    end
end

% Generate t values
for j = 0:N
    t(j+1) = j*dt;  % dt = ∆t
    % K = d*(Ev*(v-1)+4*G*vv^2*(v+1))/(Ev^2*(v-1))
end

% Main computation loop
for to = 0:N
    pin = 1*10^(-2 + to*drind);
    
    for m = 0:N
        ss = s(m+1);
        
        for n = 0:N
            tt = t(n+1);
            b(n+1) = pedelta;  % Right-hand side of the matrix equation
            
            % Calculate MeijerG function for two terms
            if pin^6/46656*(ss+tt)^6 < 10^-30
                sum1 = meijerG(1/3, [], [0,1/3,1/3,2/3,5/6], [1/6,1/2], 10^(-30)*(n+1)/(m+1));
            else
                sum1 = meijerG(1/3, [], [0,1/3,1/3,2/3,5/6], [1/6,1/2], pin^6/46656*(ss+tt)^6);
                % For incompressible layers:
                % sum1 = meijerG([1/3], [], [0,1/3,1/3,2/3], [1/6,1/2,5/6], pin^6/46656*(ss+tt)^6);
            end
            
            if pin^6/46656*(ss-tt)^6 < 10^-30
                sum2 = meijerG(1/3, [], [0,1/3,1/3,2/3,5/6], [1/6,1/2], 10^(-30)*(n+1)/(m+1));
            else
                sum2 = meijerG(1/3, [], [0,1/3,1/3,2/3,5/6], [1/6,1/2], pin^6/46656*(ss-tt)^6);
            end
            
            sum = sum1 + sum2;
            
            % Alternative method: Direct numerical integration (faster)
            % The calculation of MeijerG function in MATLAB is relatively slow.
            % To improve computing speed, consider using numerical integration directly:
            % beta1 = pin*ss + pin*tt;
            % beta2 = pin*ss - pin*tt;
            % fo = @(xi) (cos(beta1.*xi) + cos(beta2.*xi))./(xi.^3 + 1);
            % Note: pedelta should be adjusted to 2 for this method.
            % sum = integral(fo, 0, Inf);
            
            % Apply Simpson's rule coefficients
            if (n > 0) && (n < N)
                if mod(n+1, 2) == 0
                    AA(m+1, n+1) = 4*sum*dt/3;
                else
                    AA(m+1, n+1) = 2*sum*dt/3;
                end
            else
                AA(m+1, n+1) = sum*dt/3;
            end
        end  % End of n loop
    end  % End of m loop
    
    % **************************************************************
    if to == 0
        b = double(b');
    end
    
    % Solve linear system
    f = AA\b;
    F = 0;
    
    % Numerical integration using Simpson's rule
    for n = 0:N
        if (n > 0) && (n < N)
            if mod(n+1, 2) == 0
                fi(n+1) = 4*f(n+1)*dt/3;
            else
                fi(n+1) = 2*f(n+1)*dt/3;
            end
        else
            fi(n+1) = f(n+1)*dt/3;
        end
        F = F + fi(n+1);
    end
    
    Final(to+1) = F;
    Rind(to+1) = pin;
    to  % Display current iteration index
end