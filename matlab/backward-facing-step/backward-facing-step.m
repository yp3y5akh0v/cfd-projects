clear; 
clc;

nx = 150;
ny = 75;
nstepx = 15;
nstepy = 50;
lx = 4;
ly = 2;

x = linspace(0, lx, nx);
y = linspace(0, ly, ny);

xm = 0.5 * (x(1 : nx - 1) + x(2 : nx));
ym = 0.5 * (y(1 : ny - 1) + y(2 : ny));

dx = x(2) - x (1);
dy = y(2) - y (1);
dxi = 1. / dx;
dyi = 1. / dy;

nu = .001;
rho = 1;
dt = .015;

nxLaplace = nx - 1;
nyLaplace = ny - 1;

d = -2 * (dxi^2 + dyi^2);

u = zeros(nx, ny);
v = zeros(nx, ny);
us = zeros(nx, ny);
vs = zeros(nx, ny);

t = 0;
seidelItr = 100;

while t < 100

t = t + dt;

maxuv = 0;

h = y(ny) - y(nstepy) + dy;
for j = nstepy + 1 : ny - 1
    u(1, j) = - (y(j) - y(nstepy)) * ((y(j) - y(nstepy)) / h - 1);
end

v(1, nstepy + 1 : ny - 1) = 0;

u(:, ny) = 0;
v(:, ny) = 0;

u(1 : nstepx, nstepy) = 0;
v(1 : nstepx, nstepy) = 0;

u(nstepx, 1 : nstepy) = 0;
v(nstepx, 1 : nstepy) = 0;

u(nstepx : nx, 1) = 0;
v(nstepx : nx, 1) = 0;

u(nx, :) = 0;
v(nx, :) = 0;

u(1 : nstepx - 1, 1 : nstepy - 1) = 0;
v(1 : nstepx - 1, 1 : nstepy - 1) = 0;

for i = 1 : nx - 1 : nx
    for j = 1 : ny
        curuv = max(u(i, j)^2 + v(i, j)^2);
        if curuv > maxuv
            maxuv = curuv;
        end           
    end
end

for j = 1 : ny - 1 : ny
    for i = 1 : nx
        curuv = max(u(i, j)^2 + v(i, j)^2);
        if curuv > maxuv
            maxuv = curuv;
        end           
    end
end

for i = 2 : nx - 1  
    if i <= nstepx 
        j = nstepy + 1;
    else
        j = 2;
    end    
    for j = j : ny - 1      
        v_here = 0.25 * (v(i - 1, j) + v(i - 1, j + 1) + v(i, j) + v(i, j + 1));
        u_here = 0.25 * (u(i, j - 1) + u(i, j) + u(i + 1, j - 1) + u(i + 1, j));
        us(i, j) = u(i, j) + dt * ...
            (nu * (u(i - 1, j) - 2 * u(i, j) + u(i + 1, j)) * dxi^2 ...
            +nu * (u(i, j - 1) - 2 * u(i, j) + u(i, j + 1)) * dyi^2 ...
            -u(i, j) * (u(i + 1, j) - u(i - 1, j)) * 0.5 * dxi ...
            -v_here * (u(i, j + 1) - u(i, j - 1)) * 0.5 * dyi);
        vs(i, j) = v(i, j) + dt * ...
            (nu * (v(i - 1, j) - 2 * v(i, j) + v(i + 1, j)) * dxi^2 ...
            +nu * (v(i, j - 1) - 2 * v(i, j) + v(i, j + 1)) * dyi^2 ...
            -u_here  * (v(i + 1, j) - v(i - 1, j)) * 0.5 * dxi ...
            -v(i, j) * (v(i, j + 1) - v(i, j - 1)) * 0.5 * dyi);       
    end
end

p = zeros(nxLaplace, nyLaplace);
pv = zeros(nxLaplace, nyLaplace);
r = zeros(nxLaplace, nyLaplace);

for i = 1 : nxLaplace
    if i <= nstepx 
        j = nstepy + 1;
    else
        j = 2;
    end    
    for j = j : nyLaplace
        r(i, j) = rho / dt * ((us(i + 1, j) - us(i, j)) * dxi ...
                           +(vs(i, j + 1) - vs(i, j)) * dyi);        
    end
end

for k = 1 : seidelItr
    
    p(1, nstepy + 1 : nyLaplace - 1) = p(2, nstepy + 1 : nyLaplace - 1);
    p(:, nyLaplace) = p(:, nyLaplace - 1);
    p(1 : nstepx - 1, nstepy) = p(1 : nstepx - 1, nstepy + 1);
    p(nstepx, 2 : nstepy - 1) = p(nstepx + 1, 2 : nstepy - 1);
    p(nstepx + 1 : nxLaplace, 1) = p(nstepx + 1 : nxLaplace, 2);
    p(nxLaplace, :) = p(nxLaplace - 1, :);
    p(nstepx, nstepy) = 0.5 * (p(nstepx - 1, nstepy) + p(nstepx, nstepy - 1));    
    p(nstepx, 1) = 0.5 * (p(nstepx, 2) + p(nstepx + 1, 1));
    p(1 : nstepx - 1, 1 : nstepy - 1) = 0;   
    pv = p;
    
    for i = 2 : nxLaplace - 1
        if i <= nstepx 
            j = nstepy + 1;
        else
            j = 2;
        end
        for j = j : nyLaplace - 1        
            ki1 = 0;
            if i - 1 > 0
                ki1 = p(i - 1, j);
            end

            ki2 = 0;
            if i + 1 <= nxLaplace
                ki2 = pv(i + 1, j);
            end  
            
            kj1 = 0;
            if j - 1 > 0
                kj1 = p(i, j - 1);
            end

            kj2 = 0;
            if j + 1 <= nyLaplace
                kj2 = pv(i, j + 1);
            end              
            
            p(i, j) = (r(i, j) - ki1 * dxi^2 - ki2 * dxi^2 - kj1 * dyi^2 - kj2 * dyi^2) / d;
        end        
    end
    
end

for i = 2 : nx - 1
    if i <= nstepx 
        j = nstepy + 1;
    else
        j = 2;
    end    
    for j = j : ny - 1
        u(i, j) = us(i, j) - dt / rho * (p(i, j) - p(i - 1, j)) * dxi; 
        v(i, j) = vs(i, j) - dt / rho * (p(i, j) - p(i, j - 1)) * dyi;
        curuv = max(u(i, j)^2 + v(i, j)^2);
        if curuv > maxuv
            maxuv = curuv;
        end        
    end        
end

dt = min(0.25 * max(dx, dy)^2 / nu, 2 * nu / maxuv);

imagesc(xm, ym, p');
colorbar;
hold on;
quiver(x(1 : 2 : nx), y(1 : 2 : ny), u(1 : 2 : nx, 1 : 2 : ny)', v(1 : 2 : nx, 1 : 2 : ny)', 0, 'color', [0 0 0]);
axis equal;
set(gca,'YDir','normal');
set(gcf, 'Position', [480, 135, 1024, 780]);
axis([0 lx 0 ly]);
hold off;
pause(dt);

end