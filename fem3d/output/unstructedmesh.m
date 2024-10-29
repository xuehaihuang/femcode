% Define node coordinates (with boundary node information in last column)
nodes = [
    0 0 0; 1 0 0; 1 1 0; 0 1 0;  % Dirichlet nodes
    0 0 1; 1 0 1; 1 1 1; 0 1 1;  % Dirichlet nodes
    0.5 0.5 0; 0.5 0 0; 0.5 0.5 0.5;  % Inner node
    0.5 0 0.5; 0 0 0.5; 0 0.5 0.5;  % Dirichlet nodes
    0 0.5 0; 1 0.5 0.5; 1 0.5 0;  % Dirichlet nodes
    1 0 0.5; 1 1 0.5; 0.5 1 0;  % Dirichlet nodes
    0.5 1 0.5; 0 1 0.5; 0.5 0.5 1;  % Dirichlet nodes
    0.5 0 1; 0 0.5 1; 1 0.5 1;  % Dirichlet-Neumann nodes
    0.5 1 1  % Dirichlet-Neumann node
];

% Define the scale of random perturbation
scale = 0.2;  % Adjust this value for more or less perturbation

% Apply random perturbation to the nodes
random_displacement = scale * (rand(size(nodes)) - 0.5)  % Random values between -scale/2 and scale/2
II= abs(nodes-0.5)<0.1;
nodes(II) = nodes(II) + random_displacement(II);  % New node positions

nodes= [         0         0         0
    1.0000         0         0
    1.0000    1.0000         0
         0    1.0000         0
         0         0    1.0000
    1.0000         0    1.0000
    1.0000    1.0000    1.0000
         0    1.0000    1.0000
    0.5029    0.5978         0
    0.5769         0         0
    0.5176    0.5731    0.5725
    0.4310         0    0.4969
         0         0    0.5690
         0    0.5055    0.4419
         0    0.4959         0
    1.0000    0.5603    0.5260
    1.0000    0.4456         0
    1.0000         0    0.5229
    1.0000    1.0000    0.4725
    0.4180    1.0000         0
    0.4223    1.0000    0.4979
         0    1.0000    0.4385
    0.5357    0.5172    1.0000
    0.4990         0    1.0000
         0    0.5333    1.0000
    1.0000    0.4167    1.0000
    0.4295    1.0000    1.0000]

% Define element connectivity (4 nodes per element, tetrahedral)
elements = [
    7 11 19 16; 11 1 9 10; 19 9 3 17; 16 10 17 2;
    19 9 17 16; 19 9 16 11; 16 10 11 9; 16 10 9 17;
    7 11 16 26; 11 1 10 12; 16 10 2 18; 26 12 18 6;
    16 10 18 12; 16 10 12 11; 26 12 11 16; 26 12 16 18;
    7 11 26 23; 11 1 12 13; 26 12 6 24; 23 13 24 5;
    26 12 24 23; 26 12 23 11; 23 13 11 12; 23 13 12 24;
    7 11 23 27; 11 1 13 14; 23 13 5 25; 27 14 25 8;
    23 13 25 14; 23 13 14 11; 27 14 11 23; 27 14 23 25;
    7 11 27 21; 11 1 14 15; 27 14 8 22; 21 15 22 4;
    27 14 22 21; 27 14 21 11; 21 15 11 14; 21 15 14 22;
    7 11 21 19; 11 1 15 9; 21 15 4 20; 19 9 20 3;
    21 15 20 9; 21 15 9 11; 19 9 11 21; 19 9 21 20
    % Add more elements as needed...
];

% Extract triangular faces from each tetrahedral element
faces = [];
for i = 1:size(elements, 1)
    elem = elements(i, :);
    % Extract the 4 triangular faces of the tetrahedron
    faces = [faces;
        elem([1 2 3]);  % Face 1
        elem([1 2 4]);  % Face 2
        elem([1 3 4]);  % Face 3
        elem([2 3 4])   % Face 4
    ];
end

% Sort the nodes in each face and count the occurrences of each face
sorted_faces = sort(faces, 2);
[unique_faces, ~, idx] = unique(sorted_faces, 'rows', 'stable');  % 'stable' keeps the original order
face_counts = accumarray(idx, 1);  % Count occurrences of each unique face

% Determine which faces are on the boundary (appears only once)
boundary_faces = unique_faces(face_counts == 1, :);
internal_faces = unique_faces(face_counts > 1, :);

% Display the generated face information with boundary indicator
% disp('Generated face information with boundary indication:');
% for i = 1:size(unique_faces, 1)
%     if face_counts(i) == 1
%         disp(['Face ', num2str(unique_faces(i, :)), ' is on the boundary']);
%     else
%         disp(['Face ', num2str(unique_faces(i, :)), ' is internal']);
%     end
% end

% Plot the mesh with edges of boundary and internal faces
figure;
hold on;

% Plot edges of each face based on their classification
for i = 1:size(unique_faces, 1)
    f = unique_faces(i, :);
    % Define the vertices for the face
    vertices = nodes(f, :);
    
    % Extract the edges of the face
    edges = [f(1) f(2); f(2) f(3); f(3) f(1)];
    
    % Check if it's a boundary face
    if face_counts(i) == 1
        % Plot boundary face with solid line
        for j = 1:size(edges, 1)
            plot3(nodes(edges(j, :), 1), nodes(edges(j, :), 2), nodes(edges(j, :), 3), 'k-', 'LineWidth', 1.5);
        end
    else
        % Plot internal face with dashed line
        for j = 1:size(edges, 1)
            plot3(nodes(edges(j, :), 1), nodes(edges(j, :), 2), nodes(edges(j, :), 3), 'k--', 'LineWidth', 1.5);
        end
    end
end

% Mark nodes with indices
% for i = 1:size(nodes, 1)
%     text(nodes(i, 1), nodes(i, 2), nodes(i, 3), num2str(i), 'Color', 'b', 'FontSize', 8);
% end

% % xlabel('X');
% % ylabel('Y');
% % zlabel('Z');
% % title('Mesh Visualization with Boundary and Internal Faces');
% % grid on;
hold off;

az = 15; el = 30;
view(az,el)
axis square
axis off
