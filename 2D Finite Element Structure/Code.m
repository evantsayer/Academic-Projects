clear variables;
%% Pre Processing %%
file = "R";
% Read Nodes %
nodes = readmatrix("nodes" + file + ".txt");
nnodes = length(nodes(:, 1));       % Number of nodes
node = nodes(:, 2:3);               % Node matrix

% Read Elements %
elements = readmatrix("elements" + file + ".txt");
nelements = length(elements(:, 1)); % Number of elements
element = elements(:, 2:4);         % Element matrix
E = 1.0;                         % Element Young's modulus
poisson = 0.3;

% Read Forces %
forces = readmatrix("forces" + file + ".txt");
nforces = length(forces(:, 1));
fbcnd = forces(:, 1);               % Force nodes matrix
fbcdof = forces(:, 2);              % Force degrees of freedom matrix
fbcval = forces(:, 3);              % Force values matrix

% Read Displacements %
displacements = readmatrix("displacements" + file + ".txt");
ndisplacements = length(displacements(:, 1));   % Number of displacements
dbcnd = displacements(:, 1);        % Displacement nodes matrix
dbcdof = displacements(:, 2);       % Displacement degrees of freedom matrix
dbcval = displacements(:, 3);       % Displacement values matrix

% Set Constants %
r = 1;                              % Radius of Hole in structure
ndim = 2;                           % Number of dimensions
nnds = 3;                           % Number of nodes per element
ndofs = ndim * nnodes;              % Number of degrees of freedom

% Construct gcon Matrix %
gcon = zeros(nnodes, ndim);
for i = 1:nnodes
    for j = 1:ndim
        gcon(i, j) = 2 * (i - 1) + j;
    end
end

% Expand gcon Matrix
for i = 1:ndisplacements
    dofnum = gcon(dbcnd(i), dbcdof(i));
    for j = 1:nnodes
        for k = 1:ndim
            if (gcon(j, k) > dofnum)
                gcon(j, k) = gcon(j, k) - 1;
            end
        end
    end
    gcon(dbcnd(i), dbcdof(i)) = nnodes * ndim;
    ndofs = ndofs - 1;
end

%% Assembly %%
% Construct U Matrix %
U = zeros(nnodes, ndim);
for i = 1:ndisplacements
    U(dbcnd(i), dbcdof(i)) = dbcval(i);
end

% Construct Fglobal Matrix %
Fglobal = zeros(ndofs, 1);
for i = 1:nforces
    gdof = gcon(fbcnd(i), fbcdof(i));
    Fglobal(gdof) = Fglobal(gdof) + fbcval(i);
end

% Reduce Fglobal %

% Construct Kglobal Matrix %
Kglobal = zeros(ndofs, ndofs);
for iele = 1:nelements
    Kele = calculateKele(iele, E, poisson, node, element);
    for i = 1:nnds
        for j = 1:ndim
            local1 = ndim * (i - 1) + j;
            global1 = gcon(element(iele, i), j);
            if (global1 <= ndofs)
                for k = 1:nnds
                    for l = 1:ndim
                        local2 = ndim * (k - 1) + l;
                        global2 = gcon(element(iele, k), l);
                        if (global2 > ndofs)
                            Fglobal(global1) = Fglobal(global1) - Kele(local1, local2) * U(element(iele, k), l);
                        else
                            Kglobal(global1, global2) = Kglobal(global1, global2) + Kele(local1, local2);
                        end
                    end
                end
            end
        end
    end
end

%% Post-Processing %%
% Solve for Displacements %
Ureduced = linsolve(Kglobal, Fglobal);
disp("Displacements: ")
disp(Ureduced);

% Repopulate U with calculated displacements %
for i = 1:nnodes
    for j = 1:ndim
        dof = gcon(i,j);
        if (dof<=ndofs)
            U(i,j) = Ureduced(dof);
        end
    end
end

% Find Strain and Stress %
strain = zeros(nelements, 3);
stress = zeros(nelements, 3);
for i = 1:nelements
    node1 = element(i, 1);
    node2 = element(i, 2);
    node3 = element(i, 3);

    x1 = node(node1, 1);
    y1 = node(node1, 2);
    x2 = node(node2, 1);
    y2 = node(node2, 2);
    x3 = node(node3, 1);
    y3 = node(node3, 2);

    A = .5*(x2*y3-x3*y2+x3*y1-x1*y3+x1*y2-x2*y1);
 
    b1 = (.5/A)*(y2-y3);
    c1 = (.5/A)*(x3-x2);
    b2 = (.5/A)*(y3-y1);
    c2 = (.5/A)*(x1-x3);
    b3 = (.5/A)*(y1-y2);
    c3 = (.5/A)*(x2-x1);

    B = [
        b1 0 b2 0 b3 0;
        0 c1 0 c2 0 c3;
        c1 b1 c2 b2 c3 b3;
        ];


    C = [
        E/(1-poisson.^2) poisson*E/(1-poisson.^2) 0;
        poisson*E/(1-poisson.^2) E/(1-poisson.^2) 0;
        0 0 E/(2*(1+poisson));
        ];
    
    Uele = [
        U(node1, 1);
        U(node1, 2);
        U(node2, 1);
        U(node2, 2);
        U(node3, 1);
        U(node3, 2);
        ];

    strainele = B*Uele;
    strain(i, 1) = strainele(1);
    strain(i, 2) = strainele(2);
    strain(i, 3) = strainele(3);

    stressele = C*strainele;
    stress(i, 1) = stressele(1);
    stress(i, 2) = stressele(2);
    stress(i, 3) = stressele(3);
end
disp("Strains: ")
disp(strain);
disp("Stresses: ")
disp(stress);

xnode = node(:, 2) == 0;                                                % Nodes on the x axis (y = 0)
ynode = node(:, 1) == 0;                                                % Nodes on the y axis (x = 0)

xelement = sum(xnode(element), 2) >= 2;                                 % Elements on the x axis (y = 0)
yelement = sum(ynode(element), 2) >= 2;                                 % Elements on the y axis (x = 0)

nxelements = size(element(xelement, :));                                % Number of elements on the x axis (y = 0)
nyelements = size(element(yelement, :));                                % Number of elements on the y axis (x = 0)

xpos = mean(reshape(node(element(xelement, :), 1), nxelements), 2);     % midpoint of x values of elements on x axis
ypos = mean(reshape(node(element(yelement, :), 2), nyelements), 2);     % midpoint of y values of elements on y axis

xstressonx = stress(xelement, 1);                                       % x direction stress on the x axis elements
ystressonx = stress(xelement, 2);                                       % y direction stress on the x axis elements
xstressony = stress(yelement, 1);                                       % x direction stress on the y axis elements
ystressony = stress(yelement, 2);                                       % y direction stress on the y axis elements

%% Plotting %%
% Mesh %
figure(1)
hold on
title("Original Structure")
triplot(element, node(:, 1), node(:, 2), 'k')
legend('Original', 'Location', 'northeast')
xlabel('\bfx')
ylabel('\bfy')

figure(2)
hold on
title("Deformed Structure")
triplot(element, node(:, 1) + 0.1.*U(:, 1), node(:, 2) + 0.1.*U(:, 2), 'b')
legend('Deformed', 'Location', 'northeast')
xlabel('\bfx')
ylabel('\bfy')

figure(3)
hold on
title("Overlayed Deformed Structure on Original Structure")
triplot(element, node(:, 1), node(:, 2), 'k')
triplot(element, node(:, 1) + 0.1.*U(:, 1), node(:, 2) + 0.1.*U(:, 2), 'b')
legend('Original', 'Deformed', 'Location', 'northeast')
xlabel('\bfx')
ylabel('\bfy')

% Edge Stresses %
figure(4)
tiledlayout(2,2);
nexttile
interceptxx = polyval(polyfit(xpos(1:2), xstressonx(1:2), 1), r);
hold on
title("Stress X on X Axis Elements")
scatter(r, interceptxx, [], "r", "filled")
scatter(xpos, xstressonx, [], "blue", "filled")
legend('Extrapolated', 'Calculated', 'Location', 'northeast')
xlabel('\bfx')
ylabel('\bf\sigma_{xx}')

nexttile
interceptyx = polyval(polyfit(xpos(1:2), ystressonx(1:2), 1), r);
hold on
title("Stress Y on X Axis Elements")
scatter(r, interceptyx, [], "r", "filled")
scatter(xpos, ystressonx, [], "blue", "filled")
legend('Extrapolated', 'Calculated', 'Location', 'northeast')
xlabel('\bfx')
ylabel('\bf\sigma_{yy}')

nexttile
interceptxy = polyval(polyfit(ypos(1:2), xstressony(1:2), 1), r);
hold on
title("Stress X on Y Axis Elements")
scatter(r, interceptxy, [], "r", "filled")
scatter(ypos, xstressony, [], "blue", "filled")
legend('Extrapolated', 'Calculated', 'Location', 'southeast')
xlabel('\bfy')
ylabel('\bf\sigma_{xx}')

nexttile
interceptyy = polyval(polyfit(ypos(1:2), ystressony(1:2), 1), r);
hold on
title("Stress Y on Y Axis Elements")
scatter(r, interceptyy, [], "r", "filled")
scatter(ypos, ystressony, [], "blue", "filled")
legend('Extrapolated', 'Calculated', 'Location', 'southeast')
xlabel('\bfy')
ylabel('\bf\sigma_{yy}')


%% Functions %%
% Calculate Kele Function
function Kele = calculateKele(iele, E, poisson, node, element)
    x1 = node(element(iele, 1), 1);
    y1 = node(element(iele, 1), 2);
    x2 = node(element(iele, 2), 1);
    y2 = node(element(iele, 2), 2);
    x3 = node(element(iele, 3), 1);
    y3 = node(element(iele, 3), 2);

    A = .5*(x2*y3-x3*y2+x3*y1-x1*y3+x1*y2-x2*y1);

    b1 = (.5/A)*(y2-y3);
    c1 = (.5/A)*(x3-x2);
    b2 = (.5/A)*(y3-y1);
    c2 = (.5/A)*(x1-x3);
    b3 = (.5/A)*(y1-y2);
    c3 = (.5/A)*(x2-x1);

    B = [
        b1 0 b2 0 b3 0;
        0 c1 0 c2 0 c3;
        c1 b1 c2 b2 c3 b3;
        ];

    C = [
        E/(1-poisson.^2) poisson*E/(1-poisson.^2) 0;
        poisson*E/(1-poisson.^2) E/(1-poisson.^2) 0;
        0 0 E/(2*(1+poisson));
        ];
    Kele = A * B' * C * B;
end