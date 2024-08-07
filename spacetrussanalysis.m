clc;clear all;close all;
%% % Input Data
NN = xlsread('Truss_input.xls','Trussinput','F3');%Number of Nodes
Coor = xlsread('Truss_input.xls','Trussinput','D9:F100');%Coordinate
for j = 1:NN 
    i = j;
    X(i) = Coor(j,1);Y(i) = Coor(j,2);Z(i) = Coor(j,3);
end
NE = xlsread('Truss_input.xls','Trussinput','F4');%Number of elements
Pro = xlsread('Truss_input.xls','Trussinput','H9:L100');%Connectivity
for j = 1:NE
    i = Pro(j,1);E(i) = Pro(j,2);A(i) = Pro(j,3);
    Dir(i,1) = i;Dir(i,2) = Pro(j,4);Dir(i,3) = Pro(j,5);
end
NC = xlsread('Truss_input.xls','Trussinput','R3');%No. of constraints
Cons = xlsread('Truss_input.xls','Trussinput','P9:Q141');
NF = xlsread('Truss_input.xls','Trussinput','Y3');
Forc = xlsread('Truss_input.xls','Trussinput','V9:Y100');
elementNodes = xlsread('Truss_input.xls','Trussinput','K9:L100');
%% % Define Global Stiffness Matrix Element for each element
for i = 1:NE
    L(i) = sqrt((X(Dir(i,3)) - X(Dir(i,2)))^2 + (Y(Dir(i,3)) - Y(Dir(i,2)))^2 + (Z(Dir(i,3)) - Z(Dir(i,2)))^2);
    Cx(i) = (X(Dir(i,3)) - X(Dir(i,2))) / L(i);
    Cy(i) = (Y(Dir(i,3)) - Y(Dir(i,2))) / L(i);
    Cz(i) = (Z(Dir(i,3)) - Z(Dir(i,2))) / L(i);
end
%% plot
figure;
hold on;
for i = 1:NE% Plot elements
    node1 = elementNodes(i, 1);node2 = elementNodes(i, 2);
    plot3([Coor(node1, 1), Coor(node2, 1)], [Coor(node1, 2), Coor(node2, 2)], [Coor(node1, 3), Coor(node2, 3)], 'b');
end
scatter3(Coor(:, 1), Coor(:, 2), Coor(:, 3), 'r', 'filled');
for i = 1:NN% Display node coordinates
    text(Coor(i, 1), Coor(i, 2), Coor(i, 3), sprintf('(%g, %g, %g)', Coor(i, 1), Coor(i, 2), Coor(i, 3)), 'FontSize', 8);
end
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Nodes,Elements');grid on; axis equal;
%% % Initialize the stiffness matrix

K = zeros(3*NN,3*NN);
for i = 1:NE % Loop to calculate k1, k2, and k3 matrices
    k_i = (A(i) * E(i) / L(i)) * [Cx(i)^2, Cx(i)*Cy(i), Cx(i)*Cz(i), -Cx(i)^2, -Cx(i)*Cy(i), -Cx(i)*Cz(i);
                                   Cx(i)*Cy(i), Cy(i)^2, Cy(i)*Cz(i), -Cx(i)*Cy(i), -Cy(i)^2, -Cy(i)*Cz(i);
                                   Cx(i)*Cz(i), Cy(i)*Cz(i), Cz(i)^2, -Cx(i)*Cz(i), -Cy(i)*Cz(i), -Cz(i)^2;
                                  -Cx(i)^2, -Cx(i)*Cy(i), -Cx(i)*Cz(i), Cx(i)^2, Cx(i)*Cy(i), Cx(i)*Cz(i);
                                  -Cx(i)*Cy(i), -Cy(i)^2, -Cy(i)*Cz(i), Cx(i)*Cy(i), Cy(i)^2, Cy(i)*Cz(i);
                                  -Cx(i)*Cz(i), -Cy(i)*Cz(i), -Cz(i)^2, Cx(i)*Cz(i), Cy(i)*Cz(i), Cz(i)^2];
    indice = elementNodes(i, :);% Element degrees of freedom (Dof)
    elementDof = [3*indice(1)-2, 3*indice(1)-1, 3*indice(1), 3*indice(2)-2, 3*indice(2)-1, 3*indice(2)];
    K(elementDof, elementDof) = K(elementDof, elementDof) + k_i;% Assemble the global stiffness matrix
end
disp('Global Stiffness(kN/m):');disp(K);

%% % Calculate the bandwidth of the stiffness matrix K
HBW = halfBandwidthMatrix(K, 1);
disp('Half bandwidth matrix:');
disp(HBW);
%% % Definition of Primary Nodal Forces
F = zeros(3*NN,1);
for i = 1:NF
    f = 3*Forc(i,1);
    F(f-2,1) = Forc(i,2);
    F(f-1,1) = Forc(i,3);
    F(f,1) = Forc(i,4);
end
S = K;% Elimination of rows and columns of K-matrix with respect to Supports
for i = 1:NC
    r = 3*Cons(i,1);
    if Cons(i,2) == 0
        S(r-2,:) = 0; S(:,r-2) = 0; S(r-2,r-2) = 1;
        S(r-1,:) = 0; S(:,r-1) = 0; S(r-1,r-1) = 1;
        S(r,:) = 0; S(:,r) = 0; S(r,r) = 1;
    elseif Cons(i,2) == 1
        S(r-2,:) = 0; S(:,r-2) = 0; S(r-2,r-2) = 1;
    elseif Cons(i,2) == 2
        S(r-1,:) = 0; S(:,r-1) = 0; S(r-1,r-1) = 1;
    elseif Cons(i,2) == 3
        S(r,:) = 0; S(:,r) = 0; S(r,r) = 1;
    end
end

%% Gauss elimination method(Remove comment to use Gauss elimination)
GM = rref([S F]);
d = GM(:, end);
disp(d);
disp('Displacements (m):');% Display displacements
for i = 1:NN
     fprintf('dx%g',i); fprintf('=%g\n',d(3*i-2))
     fprintf('dy%g',i); fprintf('=%g\n',d(3*i-1))
     fprintf('dz%g',i); fprintf('=%g\n',d(3*i))
end
disp('-------------------------------')
disp('Support Reactions:(kN)')
%% % Calculation of Nodal Reaction, stress and strain
W = K * d;
for i = 1:NN
    fprintf('Wx%g',i); fprintf('=%g\n',W(3*i-2))
    fprintf('Wy%g',i); fprintf('=%g\n',W(3*i-1))
    fprintf('Wz%g',i); fprintf('=%g\n',W(3*i))
end
disp('-------------------------------')
disp('Stress (kPa):')
ff = zeros(NE, 6);
for e = 1:NE
    indice = elementNodes(e, :);
    elementDof = [3 * indice(1) - 2, 3 * indice(1) - 1, 3 * indice(1), 3 * indice(2) - 2, 3 * indice(2) - 1, 3 * indice(2)]; 
    x1 = Coor(indice(1), 1); y1 = Coor(indice(1), 2); z1 = Coor(indice(1), 3);
    x2 = Coor(indice(2), 1); y2 = Coor(indice(2), 2); z2 = Coor(indice(2), 3);
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
    Cx = (x2 - x1) / L; Cy = (y2 - y1) / L; Cz = (z2 - z1) / L;
    u = d(elementDof);
    memberstress(e) = E(e) / L * [-Cx, -Cy, -Cz, Cx, Cy, Cz] * u;
    fprintf('%3d %12.8f\n', e, memberstress(e));
end
disp('-------------------------------')
disp('Elements Strain:')
for i = 1:NE
    strain(i)=memberstress(i)/E(i);
    fprintf('Strain %g:\t', i); 
    fprintf('%g\n', strain(i));
end
%Write the matrices to Excel
xlswrite('Trussresult2.xlsx', K, 'Stiffness_Matrix');
xlswrite('Trussresult2.xlsx', d, 'Displacements');
xlswrite('Trussresult2.xlsx', W,'Support_Reactions');
xlswrite('Trussresult2.xlsx', memberstress', 'Member_Stresses');
xlswrite('Trussresult2.xlsx', strain', 'Member_Strain');
xlswrite('Trussresult2.xlsx', HBW, 'Half_band_width_Matrix');
disp('Results have been successfully exported to Excel.');


%% Function for Half band width matrix
function HBW = halfBandwidthMatrix(K, half_bandwidth)
    % Get the size of matrix A
    [m, n] = size(K);

    % Initialize a zero matrix for the half-bandwidth matrix
    HBW = zeros(m, half_bandwidth*2+1);

    % Loop through the elements of A to form the half-bandwidth matrix
    for i = 1:m
        for j = max(1, i - half_bandwidth):min(n, i + half_bandwidth)
            HBW(i, j - (i - half_bandwidth) + 1) = K(i, j);
        end
    end
end
