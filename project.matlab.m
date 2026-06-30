%% this is the project code by group 3
% please note that we have used the approximation of bessel function for
% the calculation of c(x) as we are not allowed to use the exact and better
% in built subroutines available in matlab like besseli(),etc. 
% please also note that we have also submitted the code in which we use the
% inbuilt matlab subroutine for analytical solution i.e. use the besseli
% funns
%% members-----
% 2024chb1066- dhruv kaushal
% 2024chb1062- cheemakurthi purandhar sai
% 2024chb1067- ganedi durga tejaswi
% 2024chb1072- minaksh sharma
% 2024chb1086- sunitha nayak
%% we have presented the complete code as a script file instead of writing it in different discrete functions, we have added comments for explaining the code 
clc; clear; close all;
fprintf('1.lets start the work...\n');
%% given data constants----------------
L = 5;                  
dx = 0.05;              
D0 = 0.84;         
Dn = 1.09;           
k_rate = 0.1; % rate constant          
CA0 = 100;         
CAn = 100;           
%% calculating the other variables 
N = round(L/dx) + 1; % //calculating the no. of points used in numerical calculations
x = linspace(0, L, N)'; % doing transpose of linspace
slope = (Dn-D0)/ L; %calculation of slope from diffusivities
D_x = slope*x + D0;% //
dDdx = slope; %//         
%% Finite Difference Method (using guass elimination with partial pivoting) ---------
N_in = N-2; %consider only interior points
A = zeros(N_in, N_in); % matrix of alpha,beta and gamma
b = zeros(N_in, 1);     % matrix for the rhs
for r = 1:N_in
    i=r+1; % doing only on interior points
    Di = D_x(i);  
    alpha = (Di / dx^2) - (slope/(2*dx));
    beta  = (-2*Di/dx^2) - k_rate;
    gamma = (Di/dx^2) + (slope/(2*dx));
    A(r, r) = beta;
    if r > 1; 
        A(r, r-1) = alpha; 
    end
    if r < N_in; 
        A(r, r+1) = gamma; 
    end
    if r == 1; % boundary condition for x=0
        b(r) = -alpha * CA0;
    elseif r == N_in; %boundary condition for x=n
         b(r) = -gamma * CAn; 
    end
end
% Solve Ax = b using partial pivoting guass elimination
B=[A, b]; % augmented matrix for guass elimination forward elimination
n_sys = N_in;
for k = 1:N_in-1
    [~,p] = max(abs(B(k:N_in, k))); % partial pivoting
    p_row = p + k - 1;
    if p_row ~= k; 
        temp = B(k, :); B(k, :) = B(p_row, :); B(p_row, :) = temp; 
    end
    for i = k+1:N_in
        factor = B(i, k) / B(k, k);
        B(i, k:N_in+1) = B(i, k:N_in+1) - factor * B(k, k:N_in+1); % forward elimination
    end
end
C_in = zeros(N_in, 1);  % creating concn matrix for interior pts
C_in(N_in) = B(N_in, N_in+1) / B(N_in,N_in);
for i = N_in-1:-1:1
    sum_val = sum(B(i, i+1:N_in) .* C_in(i+1:N_in)'); 
    C_in(i) = (B(i,N_in+1) - sum_val) / B(i, i);
end
C_augmented = [CA0; C_in; CAn];%creating augmented matrix of concentrations
fprintf('2.numerical data  is ready. range: [%.2f, %.2f]\n', min(C_augmented), max(C_augmented));
%% spline curves calculations
n_in = N - 2;
S = zeros(n_in, n_in); Y = zeros(n_in, 1);
for i = 1:n_in
    idx = i + 1;
    S(i, i) = 4;
    if i > 1; S(i, i-1) = 1; end;
    if i < n_in; S(i, i+1) = 1; end;
    Y(i) = (6 / dx^2) * (C_augmented(idx-1) - 2*C_augmented(idx) + C_augmented(idx+1));
end
augmentedS = [S, Y];
% guass elimination
for k=1:n_in-1
    [~, p] = max(abs(augmentedS(k:n_in, k))); %pivoting
    p_row=p+k-1;
    if p_row ~= k; 
        temp = augmentedS(k, :); 
        augmentedS(k, :) = augmentedS(p_row, :);
        augmentedS(p_row, :) = temp;
    end
    for i = k+1:n_in
        factor = augmentedS(i, k) / augmentedS(k, k);
        augmentedS(i, k:n_in+1) = augmentedS(i, k:n_in+1) - factor * augmentedS(k, k:n_in+1);
    end
end
M_in = zeros(n_in, 1); % matrix for the interior points M calculation
M_in(n_in) = augmentedS(n_in, n_in+1) / augmentedS(n_in, n_in);
for i = n_in-1:-1:1
    sum_val = sum(augmentedS(i, i+1:n_in) .* M_in(i+1:n_in)');
    M_in(i) = (augmentedS(i, n_in+1) - sum_val) / augmentedS(i, i);
end
M_augmented = [0; M_in; 0];
x_values = linspace(0, L, 200)';
C_spline = zeros(size(x_values));
for k = 1:length(x_values)
    xk = x_values(k);
    idx = floor(xk / dx)+1;
    if idx >= N; 
        idx = N - 1; 
    end
    h = dx; 
    xi = x(idx); 
    xi1 = x(idx+1);
    AA = (xi1 - xk)/h;
    BB = (xk - xi)/h;
    C_spline(k) = AA*C_augmented(idx) + BB*C_augmented(idx+1) + ((AA^3 - AA)*M_augmented(idx) + (BB^3 - BB)*M_augmented(idx+1)) * (h^2)/6;
end
%% analytical solution --------------------
% c(x)= AIo(z(x))+BKo(z(x))  , z(x)=2/beta*sqrt(k(alpha+beta*x));
% these above are derived using bessel's zeroth order equation
% A= 4.029174*(10)^(-3) , B=1.435568*(10)^(7);
% c(x)= 4.029174*(10)^(-3)Io(z(x))+ 1.435568*(10)^(7)*Ko(z(x));
% Io(x)= bessel funcn of order zero of first kind = geometric series=
% Io(z)=e^z/sqrt(2pie*z) 
% Ko(x)= bessel funn of order zeros of second kind = similarly is approx, =
% sqrt(pie/2*z) *e^(-z)

% -- INTERNAL MINI FUNCTIONS FOR APPROXIMATION (Asymptotic Series) --
% We define these handles to perform the approximation calculations
Io_approx = @(z) (exp(z) ./ sqrt(2*pi*z)) .* (1 + 1./(8*z));
Ko_approx = @(z) (sqrt(pi./(2*z)) .* exp(-z)) .* (1 - 1./(8*z));

z_func = @(v) (2/slope) * sqrt(k_rate * (slope*v + D0));

% Using the mini functions instead of built-in bessel
val_i0_0 = Io_approx(z_func(0)); val_k0_0 = Ko_approx(z_func(0));
val_i0_L = Io_approx(z_func(L)); val_k0_L = Ko_approx(z_func(L));

det_val = val_i0_0 * val_k0_L - val_k0_0 * val_i0_L;
c1 = (100 * val_k0_L - 100 * val_k0_0) / det_val;
c2 = (100 * val_i0_0 - 100 * val_i0_L) / det_val;

C_anal = c1 * Io_approx(z_func(x_values)) + c2 * Ko_approx(z_func(x_values));

% SAFETY CHECK: If C_anal has complex numbers, plot will fail. We take REAL part only.
if ~isreal(C_anal)
    fprintf('WARNING: Analytical solution has complex parts. Taking real part only.\n');
    C_anal = real(C_anal);
end
fprintf('3. All Data Calculated successfully.\n');
%% ROBUST PLOTTING ------------------------
% Figure 1: Numerical vs Analytical
figure(1); 
plot(x, C_augmented, 'bo', x_values, C_anal, 'r-'); 
title('Part (i): Numerical (Points) vs Analytical (Line)');
xlabel('Distance'); ylabel('Concentration');
legend('Numerical FDM', 'Analytical');
grid on;
% Figure 2: Spline vs Analytical
figure(2);
plot(x_values, C_spline, 'go--', 'MarkerSize', 3, 'MarkerIndices', 1:5:length(x_values)); 
hold on;
plot(x_values, C_anal, 'r-', 'LineWidth', 1.5);
title('Part (ii): Spline (circle&Dashed) vs Analytical (Solid)');
xlabel('Distance'); ylabel('Concentration');
legend('Cubic Spline', 'Analytical');
grid on;