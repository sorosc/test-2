% Pricing Asian Call Option with Locally One-Dimensional (LOD) numerical scheme

% Input
s0 = 5;
K = 5;
r = 0.1;
sigma = 0.5;
T = 1;
dt = 0.001;
ds = 0.5;
dI = 0.5;

% Locally One-Dimensional scheme
s_max = 3*s0; 
nums = round(s_max/ds);
numt = 2*round(T/dt);
I_T_max = T*s_max;
numI_max = round(I_T_max/dI);
sol = zeros(numt+1,nums+1,numI_max+1);
% 1. Determine the terminal condition   
for m = 1:numI_max+1
    for n = 1:nums+1
        sol(1,n,m) = subplus((1/T)*(m-1)*I_T_max/numI_max - K);
    end
end
% 2. Determine the far boundaries (back & right); the whole boundary
for m = 2:numt+1
    t = T*(numt-m+1)/numt; 
    I_t_max = t*s_max;
    numI = numI_max;    
    for n = 2: nums+1 
        s_t = s_max*(n-1)/nums;
        sol(m,n,numI+1) = subplus((1/T)*(I_t_max + (exp(r*(T-t))-1)*s_t/r)-K);
    end
    for h = 2:numI
        I_t = I_t_max*(h-1)/numI;
        sol(m,nums+1,h) = subplus((1/T)*(I_t+(exp(r*(T-t))-1)*s_max/r)-K);
    end     
end

for m = 2:2:numt
    t = T*(numt-m+1)/numt; 
    I_t_max = t*s_max;
    numI = numI_max;
    if numI < 1
        numI = 1;
    end
% 3. Determine the near boundaries (front & left); one layer at a time
% 3.1 Corners
    % front & right
    i = nums+1;
    j = 1;
    sol(m,i,j) = sol(m-1,i,j)+dt/2*((r*i)*(sol(m-1,i,j)-sol(m-1,i-1,j))+0.5*(sigma*i)^2*(sol(m-1,i,j)-2*sol(m-1,i-1,j)+sol(m-1,i-2,j))+i*ds/dI*(sol(m-1,i,j+1)-sol(m-1,i,j))-r*sol(m-1,i,j));
    % front & left
    i = 1;
    j = 1;
    sol(m,i,j) = sol(m-1,i,j)+dt/2*((r*i)*(sol(m-1,i+1,j)-sol(m-1,i,j))+0.5*(sigma*i)^2*(sol(m-1,i+2,j)-2*sol(m-1,i+1,j)+sol(m-1,i,j))+i*ds/dI*(sol(m-1,i,j+1)-sol(m-1,i,j))-r*sol(m-1,i,j));
    % back & left
    i = 1;
    j = numI+1;
    sol(m,i,j) = sol(m-1,i,j)+dt/2*((r*i)*(sol(m-1,i+1,j)-sol(m-1,i,j))+0.5*(sigma*i)^2*(sol(m-1,i+2,j)-2*sol(m-1,i+1,j)+sol(m-1,i,j))+i*ds/dI*(sol(m-1,i,j)-sol(m-1,i,j-1))-r*sol(m-1,i,j));
% 3.2 Interior boundary
    % S-direction
    for n = 2:nums
        i = n;
        j = 1;
        sol(m,i,j) = sol(m-1,i,j)+dt/2*((r*i/2)*(sol(m-1,i+1,j)-sol(m-1,i-1,j))+0.5*(sigma*i)^2*(sol(m-1,i+1,j)-2*sol(m-1,i,j)+sol(m-1,i-1,j))+i*ds/dI*(sol(m-1,i,j+1)-sol(m-1,i,j))-r*sol(m-1,i,j));
    end
    % I-direction
    for h = 2:numI
        i = 1;
        j = h;
        sol(m,i,j) = sol(m-1,i,j)+dt/2*((r*i)*(sol(m-1,i+1,j)-sol(m-1,i,j))+0.5*(sigma*i)^2*(sol(m-1,i+2,j)-2*sol(m-1,i+1,j)+sol(m-1,i,j))+i*ds/(dI*2)*(sol(m-1,i,j+1)-sol(m-1,i,j-1))-r*sol(m-1,i,j));
    end
% 4. Calculate the interior of an intermediate layer (n+0.5); fully implicit; expanding in the s-direction
    for h = 2:numI
        A = zeros(nums-1,nums-1);
        b = zeros(nums-1,1);
        j = h;
        for n = 1:nums-1
            i = n;
            ai = dt*(0.5*r*i-0.5*(sigma*i)^2);
            bi = r*dt+(sigma*i)^2*dt+1;
            ci = dt*(-0.5*r*i-0.5*(sigma*i)^2);
            if n == 1
                A(n,n) = bi;
                A(n,n+1) = ci;
                b(n,1) = sol(m-1,i+1,j)-ai*sol(m,i,j);
            elseif n == nums-1
                    A(n,n-1) = ai;
                    A(n,n) = bi;
                    b(n,1) = sol(m-1,i+1,j)-ci*sol(m,i+2,j); 
                else
                    A(n,n-1) = ai;
                    A(n,n) = bi;
                    A(n,n+1) = ci;
                    b(n,1) = sol(m-1,i+1,j);
            end
        end
        [L U] = lu(A);
        y = L\b;
        x = U\y;        
        for g = 2:nums
            i = g;
            sol(m,i,j) = x(i-1,1);
        end     
    end
% 5. Determine the boundaries of layer n
m = m+1;
% 5.1 Corners
    % front & right
    i = nums+1;
    j = 1; 
    sol(m,i,j) = sol(m-1,i,j)+dt/2*((r*i)*(sol(m-1,i,j)-sol(m-1,i-1,j))+0.5*(sigma*i)^2*(sol(m-1,i,j)-2*sol(m-1,i-1,j)+sol(m-1,i-2,j))+i*ds/dI*(sol(m-1,i,j+1)-sol(m-1,i,j))-r*sol(m-1,i,j));
    % front & left
    i = 1;
    j = 1;
    sol(m,i,j) = sol(m-1,i,j)+dt/2*((r*i)*(sol(m-1,i+1,j)-sol(m-1,i,j))+0.5*(sigma*i)^2*(sol(m-1,i+2,j)-2*sol(m-1,i+1,j)+sol(m-1,i,j))+i*ds/dI*(sol(m-1,i,j+1)-sol(m-1,i,j))-r*sol(m-1,i,j));
    % back & left
    i = 1;
    j = numI+1;
    sol(m,i,j) = sol(m-1,i,j)+dt/2*((r*i)*(sol(m-1,i+1,j)-sol(m-1,i,j))+0.5*(sigma*i)^2*(sol(m-1,i+2,j)-2*sol(m-1,i+1,j)+sol(m-1,i,j))+i*ds/dI*(sol(m-1,i,j)-sol(m-1,i,j-1))-r*sol(m-1,i,j));
% 5.2 Interior boundary
    % S-direction
    for n = 2:nums
        i = n;
        j = 1;
        sol(m,i,j) = sol(m-1,i,j)+dt/2*((r*i/2)*(sol(m-1,i+1,j)-sol(m-1,i-1,j))+0.5*(sigma*i)^2*(sol(m-1,i+1,j)-2*sol(m-1,i,j)+sol(m-1,i-1,j))+i*ds/dI*(sol(m-1,i,j+1)-sol(m-1,i,j))-r*sol(m-1,i,j));
    end
    % I-direction
    for h = 2:numI
        i = 1;
        j = h;
        sol(m,i,j) = sol(m-1,i,j)+dt/2*((r*i)*(sol(m-1,i+1,j)-sol(m-1,i,j))+0.5*(sigma*i)^2*(sol(m-1,i+2,j)-2*sol(m-1,i+1,j)+sol(m-1,i,j))+i*ds/(dI*2)*(sol(m-1,i,j+1)-sol(m-1,i,j-1))-r*sol(m-1,i,j));
    end
% 6. Calculate the interior of layer n; fully implicit; expanding in the I-direction                    
    for n = 2:nums
        A = zeros(numI-1,numI-1);
        b = zeros(numI-1,1);
        i = n;
        for h = 1:numI-1
            j = h;
            di = i*dt*ds/(2*dI);
            ei = 1;
            fi = -i*ds*dt/(2*dI);
            if h == 1
                A(h,h) = ei;
                A(h,h+1) = fi;
                b(h,1) = sol(m-1,i,j+1)-di*sol(m,i,j);
            elseif h == numI-1
                    A(h,h-1) = di;
                    A(h,h) = ei;
                    b(h,1) = sol(m-1,i,j+1)-di*sol(m,i,j+2); 
                else
                    A(h,h-1) = ai;
                    A(h,h) = bi;
                    A(h,h+1) = ci;
                    b(h,1) = sol(m-1,i,j+1);
            end
        end
        [L U] = lu(A);
        y = L\b;
        x = U\y;        
        for g = 2:numI
            j = g;
            sol(m,i,j) = x(j-1,1);
        end     
    end
end

disp(sprintf('V(0,0.45,0) = %.4f',sol(numt+1,round(4.5/ds)+1,1)));             
disp(sprintf('V(0,0.5,0) = %.4f',sol(numt+1,round(5/ds)+1,1)));
disp(sprintf('V(0,0.55,0) = %.4f',sol(numt+1,round(5.5/ds)+1,1)));


