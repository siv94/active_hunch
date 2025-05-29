%%%%%%%%%%%%%% Lower bound FELA code for retaining wall with suction, Vs
%%%%%%%%%%%%%% and D varying with depth

clear all;
close all;

%%%% wall geometry
H=6;
H1=0.5*H;
H2=H-H1;
hr=H1/H2;
theta1=deg2rad(10);
theta2=atan(hr*tan(theta1));
% theta2=deg2rad(0);  % haunch angle with vertical
% theta1=atan((1/hr)*tan(theta2))
beta=deg2rad(0);

%%%%%%% mesh from ABAQUS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

filename = '6mhnch_nofondtn_2.inp';         % import mesh generated in ABAQUS as .inp file
fid = fopen(filename, 'r');

% Check if file opened successfully
if fid == -1
    error('Cannot open file: %s', filename);
end


%%% extracting node and elements data from mesh file

nodes = [];
elements = [];

% Read file line by line
while ~feof(fid)
    tline = fgetl(fid);
    
    % Check for the *Node section
    if contains(tline, '*Node')
        % Read nodes
        while true
            tline = fgetl(fid);
            if isempty(tline) || startsWith(tline, '*')
                break;
            end
            % Parse the node line
            data = sscanf(tline, '%f, %f, %f, %f');
            nodes = [nodes; data'];
        end
    end
    
    % Check for the *Element section
    if contains(tline, '*Element')
        % Read elements
        while true
            tline = fgetl(fid);
            if isempty(tline) || startsWith(tline, '*')
                break;
            end
            % Parse the element line
            data = sscanf(tline, '%f, %f, %f, %f, %f');
            elements = [elements; data'];
        end
    end
end

fclose(fid);

figure;
hold on;
axis equal;

%%%% TO separate nodes for soil and wall
split_index = find(nodes(2:end, 1) < nodes(1:end-1, 1), 1, 'first') + 1;

nodes_part1 = nodes(1:split_index-1, :);    % Nodes for first part (wall or soil, check in .inp)
nodes_part2 = nodes(split_index:end, :);    % Nodes for second part (wall or soil, check in .inp)

split_index1 = find(elements(2:end, 1) < elements(1:end-1, 1), 1, 'first') + 1;

T_soil=elements(1:split_index1-1, :);
T_wall=elements(split_index1:end, :);
T1=T_wall(:,2:4);       % element connectvity matrix for wall nd soil (check which part is wall an which is soil in.inp)
T2=T_soil(:,2:4);

x_wall=nodes_part2(:,2);
y_wall=nodes_part2(:,3);
x_soil=nodes_part1(:,2);
y_soil=nodes_part1(:,3);

%%%% extend of domain %%%
dom_W=max(nodes_part1(:,2));

triplot(T1,x_wall,y_wall,"-r")
triplot(T2,x_soil,y_soil,"-k")

centroids_s = zeros(size(T2, 1), 2);      % check T1 or T2 is soil, change accordingly
centroids_w = zeros(size(T1, 1), 2);

for i = 1:size(T2, 1)
       node_indices_s = T2(i, :);
       coords_s = nodes_part1(node_indices_s, 2:3);     % check which part is soil
           centroid_s = mean(coords_s, 1);
           centroids_s(i, :) = centroid_s;
end
for i = 1:size(T1, 1)
       node_indices_w = T1(i, :);
       coords_w = nodes_part2(node_indices_w, 2:3);     % check which part is soil
           centroid_w = mean(coords_w, 1);
           centroids_w(i, :) = centroid_w;
end

%%%%% identifying discontinuities

%%% dicontinuity matrix for soil
T2_renamed = zeros(size(T2));
new_node_counter = 1;
renamed_nodes_map_s = containers.Map('KeyType', 'int32', 'ValueType', 'any');

for i = 1:size(T2, 1)
    T2_renamed(i, :) = new_node_counter:new_node_counter+2;
    
    for j = 1:3
        original_node = T2(i, j);
        renamed_node_s = T2_renamed(i, j);
        renamed_nodes_map_s(renamed_node_s) = nodes_part1(original_node, 2:3);
    end
    % Increment the counter by 3 for the next set of nodes
    new_node_counter = new_node_counter + 3;
end

renamed_nodes_matrix_s = zeros(new_node_counter-1, 3);

for k = 1:new_node_counter-1
    renamed_nodes_matrix_s(k, 1) = k;  % Node number
    renamed_nodes_matrix_s(k, 2:3) = renamed_nodes_map_s(k);  % Coordinates
end

edges_map_s = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:size(T2, 1)
    tri_nodes = T2(i, :);
    coords = nodes(tri_nodes, 2:3);
    renamed_nodes_s = T2_renamed(i, :);
    tri_edges_s = [renamed_nodes_s([1, 2]);
                 renamed_nodes_s([2, 3]);
                 renamed_nodes_s([3, 1])];

    tri_edges_s = sort(tri_edges_s, 2);

    for j = 1:size(tri_edges_s, 1)
        edge_coords_s = sort(coords([j, mod(j, 3) + 1], :), 1);
        edge_key_s = mat2str(edge_coords_s(:)');

        if isKey(edges_map_s, edge_key_s)
            edges_map_s(edge_key_s) = [edges_map_s(edge_key_s); tri_edges_s(j, :)];
        else
            edges_map_s(edge_key_s) = tri_edges_s(j, :);
        end
    end
end

discontinuities_s = [];

keys_s = edges_map_s.keys;
for k = 1:length(keys_s)
    edge_nodes_s = edges_map_s(keys_s{k});
    if size(edge_nodes_s, 1) == 2
        discontinuity_row = [edge_nodes_s(1, :), edge_nodes_s(2, :)];
        discontinuities_s = [discontinuities_s; discontinuity_row];
    end
end

discont_reard_s = zeros(size(discontinuities_s));

for i = 1:size(discontinuities_s, 1)
    current_row = discontinuities_s(i, :);
    coords1 = renamed_nodes_matrix_s(current_row(1), 2:3);
    coords2 = renamed_nodes_matrix_s(current_row(2), 2:3);
    coords3 = renamed_nodes_matrix_s(current_row(3), 2:3);
    coords4 = renamed_nodes_matrix_s(current_row(4), 2:3);
    
    if isequal(coords1, coords3)
        discont_reard_s(i, :) = [current_row(1), current_row(3), current_row(2), current_row(4)];
    elseif isequal(coords1, coords4)
        discont_reard_s(i, :) = [current_row(1), current_row(4), current_row(2), current_row(3)];
    elseif isequal(coords2, coords3)
        discont_reard_s(i, :) = [current_row(2), current_row(3), current_row(1), current_row(4)];
    else
        discont_reard_s(i, :) = [current_row(2), current_row(4), current_row(1), current_row(3)];
    end
end

ddata_s=zeros(size(discont_reard_s,1),5);
ddata_s(:,1:4)=discont_reard_s;

figure;
hold on;
axis equal;

for i=1:size(ddata_s,1)
lx=renamed_nodes_matrix_s(ddata_s(i,4),2)-renamed_nodes_matrix_s(ddata_s(i,2),2);
ly=renamed_nodes_matrix_s(ddata_s(i,4),3)-renamed_nodes_matrix_s(ddata_s(i,2),3);
ddata_s(i,5)=atan(ly/lx);
xo=[renamed_nodes_matrix_s(ddata_s(i,1),2) renamed_nodes_matrix_s(ddata_s(i,4),2)];
yo=[renamed_nodes_matrix_s(ddata_s(i,1),3)  renamed_nodes_matrix_s(ddata_s(i,4),3)];
    plot(xo,yo);
end

hold off


%%%% discont matrix for wall %%%%%
T1_renamed = zeros(size(T1));
new_node_counter = 1;
renamed_nodes_map_w = containers.Map('KeyType', 'int32', 'ValueType', 'any');

% Iterate through each triangle
for i = 1:size(T1, 1)
    T1_renamed(i, :) = new_node_counter:new_node_counter+2;
    
    for j = 1:3
        original_node = T1(i, j);
        renamed_node_w = T1_renamed(i, j);
        renamed_nodes_map_w(renamed_node_w) = nodes_part2(original_node, 2:3);
    end
    new_node_counter = new_node_counter + 3;
end

renamed_nodes_matrix_w = zeros(new_node_counter-1, 3);

for k = 1:new_node_counter-1
    renamed_nodes_matrix_w(k, 1) = k;  % Node number
    renamed_nodes_matrix_w(k, 2:3) = renamed_nodes_map_w(k);  % Coordinates
end

edges_map_w = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:size(T1, 1)
    tri_nodes = T1(i, :);
    coords = nodes(tri_nodes, 2:3);
    renamed_nodes_w = T1_renamed(i, :);

    tri_edges_w = [renamed_nodes_w([1, 2]);
                 renamed_nodes_w([2, 3]);
                 renamed_nodes_w([3, 1])];

    tri_edges_w = sort(tri_edges_w, 2);

    for j = 1:size(tri_edges_w, 1)
        edge_coords_w = sort(coords([j, mod(j, 3) + 1], :), 1);
        edge_key_w = mat2str(edge_coords_w(:)');

        if isKey(edges_map_w, edge_key_w)
            edges_map_w(edge_key_w) = [edges_map_w(edge_key_w); tri_edges_w(j, :)];
        else
            edges_map_w(edge_key_w) = tri_edges_w(j, :);
        end
    end
end

discontinuities_w = [];

keys_w = edges_map_w.keys;
for k = 1:length(keys_w)
    edge_nodes_w = edges_map_w(keys_w{k});
    if size(edge_nodes_w, 1) == 2
        discontinuity_row = [edge_nodes_w(1, :), edge_nodes_w(2, :)];
        discontinuities_w = [discontinuities_w; discontinuity_row];
    end
end

discont_reard_w = zeros(size(discontinuities_w));

for i = 1:size(discontinuities_w, 1)
    current_row = discontinuities_w(i, :);
    
    coords1 = renamed_nodes_matrix_w(current_row(1), 2:3);
    coords2 = renamed_nodes_matrix_w(current_row(2), 2:3);
    coords3 = renamed_nodes_matrix_w(current_row(3), 2:3);
    coords4 = renamed_nodes_matrix_w(current_row(4), 2:3);
    
    if isequal(coords1, coords3)
        discont_reard_w(i, :) = [current_row(1), current_row(3), current_row(2), current_row(4)];
    elseif isequal(coords1, coords4)
        discont_reard_w(i, :) = [current_row(1), current_row(4), current_row(2), current_row(3)];
    elseif isequal(coords2, coords3)
        discont_reard_w(i, :) = [current_row(2), current_row(3), current_row(1), current_row(4)];
    else
        discont_reard_w(i, :) = [current_row(2), current_row(4), current_row(1), current_row(3)];
    end
end

ddata_w=zeros(size(discont_reard_w,1),5);
ddata_w(:,1:4)=discont_reard_w;

figure;
hold on;
axis equal;

for i=1:size(ddata_w,1)
lx=renamed_nodes_matrix_w(ddata_w(i,4),2)-renamed_nodes_matrix_w(ddata_w(i,2),2);
ly=renamed_nodes_matrix_w(ddata_w(i,4),3)-renamed_nodes_matrix_w(ddata_w(i,2),3);
ddata_w(i,5)=atan(ly/lx);
xo_w=[renamed_nodes_matrix_w(ddata_w(i,1),2) renamed_nodes_matrix_w(ddata_w(i,4),2)];
yo_w=[renamed_nodes_matrix_w(ddata_w(i,1),3)  renamed_nodes_matrix_w(ddata_w(i,4),3)];
    plot(xo_w,yo_w);
end

hold off

%%%% boundary matrices

extractEdges = @(T) [T(:, [1, 2]);
                     T(:, [2, 3]);
                     T(:, [3, 1])];

% Extract edges for wall and soil
wall_edges = extractEdges(T1_renamed);
soil_edges = extractEdges(T2_renamed);

% Initialize an empty array to store the shared edges
shared_edges = [];

% Loop through each edge in soil_edges and check if it matches with any wall_edge
for i = 1:size(soil_edges, 1)
    % Get the coordinates for the current soil edge
    soil_node1_coords = renamed_nodes_matrix_s(soil_edges(i, 1), 2:3);
    soil_node2_coords = renamed_nodes_matrix_s(soil_edges(i, 2), 2:3);
    
    % Loop through each wall edge to find a match
    for j = 1:size(wall_edges, 1)
        wall_node1_coords = renamed_nodes_matrix_w(wall_edges(j, 1), 2:3);
        wall_node2_coords = renamed_nodes_matrix_w(wall_edges(j, 2), 2:3);
        
        % Check if the coordinates of the soil edge match with the wall edge
        if (isequal(soil_node1_coords, wall_node1_coords) && isequal(soil_node2_coords, wall_node2_coords)) || ...
           (isequal(soil_node1_coords, wall_node2_coords) && isequal(soil_node2_coords, wall_node1_coords))
            % If a match is found, add the edge to the shared_edges array
            shared_edges = [shared_edges; soil_edges(i, :)];
            break;
        end
    end
end

shared_edge_coords = [];

for i = 1:size(shared_edges, 1)
    shrd_edg_n1=shared_edges(i, 1);
    shrd_edg_n2=shared_edges(i, 2);
    node1_coords = renamed_nodes_matrix_s(shared_edges(i, 1), 2:3);
    node2_coords = renamed_nodes_matrix_s(shared_edges(i, 2), 2:3);
    % shared_edge_coords = [shared_edge_coords; node1_coords; node2_coords];
    shared_edge_coords = [shared_edge_coords;shrd_edg_n1 node1_coords;shrd_edg_n2 node2_coords];
end
hold off;
bbw=shared_edge_coords;
% bbw=sortrows(bbw,3);
figure;
hold on;
axis equal;
for i= 1:size(bbw,1)
xb(i)=bbw(i,2);    yb(i)=bbw(i,3);
plot (xb,yb)
end

%%%%% free surface nodes

bottm_bound=[]; side_bound=[];
for i=1:size(soil_edges,1)
    edge_b=soil_edges(i,:);
nod1_y=renamed_nodes_matrix_s(edge_b(1,1),3);
nod2_y=renamed_nodes_matrix_s(edge_b(1,2),3);
nod1_x=renamed_nodes_matrix_s(edge_b(1,1),2);
nod2_x=renamed_nodes_matrix_s(edge_b(1,2),2);
nod1_x=round(nod1_x,1);
nod2_x=round(nod2_x,1);
nod1_y=round(nod1_y,1);
nod2_y=round(nod2_y,1);
if nod1_y==0 && nod2_y==0
    bottm_bound=[bottm_bound;edge_b];
end
if nod1_x==dom_W && nod2_x==dom_W
    side_bound=[side_bound;edge_b];
end
end
disc_mtx1=[discont_reard_s(:,3)  discont_reard_s(:,1);discont_reard_s(:,2)   discont_reard_s(:,4)];
disc_mtx2=[discont_reard_s(:,1)  discont_reard_s(:,3);discont_reard_s(:,4)   discont_reard_s(:,2)];
combnd_mtx=[disc_mtx1;disc_mtx2;bottm_bound;side_bound;shared_edges];

exclude_idx = false(size(soil_edges, 1), 1);
for i = 1:size(combnd_mtx, 1)
        exclude_idx = exclude_idx | ismember(soil_edges, combnd_mtx(i, :), 'rows');
end
free_bound_edgs = soil_edges(~exclude_idx, :);
bfr_s=free_bound_edgs;

%%%%% free bound of wall

bottm_bound_w=[]; top_bound_w=[];   side_bound_w=[];
for i=1:size(wall_edges,1)
    edge_w=wall_edges(i,:);
nod1_y=renamed_nodes_matrix_w(edge_w(1,1),3);
nod2_y=renamed_nodes_matrix_w(edge_w(1,2),3);
nod1_x=renamed_nodes_matrix_w(edge_w(1,1),2);
nod2_x=renamed_nodes_matrix_w(edge_w(1,2),2);
nod1_x=round(nod1_x,1);
nod2_x=round(nod2_x,1);
nod1_y=round(nod1_y,1);
nod2_y=round(nod2_y,1);
if nod1_y==0 && nod2_y==0
    bottm_bound_w=[bottm_bound_w;edge_w];
end
if nod1_y==H && nod2_y==H
    top_bound_w=[top_bound_w;edge_w];
end
if nod1_x==0 && nod2_x==0
    side_bound_w=[side_bound_w;edge_w];
end
end

bfr_w=top_bound_w;
bfr_wside=side_bound_w;


%%%%% Main program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=size(renamed_nodes_matrix_s,1)+size(renamed_nodes_matrix_w,1);
E=round(N/3);
gammaD=20;
G=2.7;
gamma_w=25;
c=0.00001;
cw=4210;
phi_w=deg2rad(40.7);
phi=deg2rad(30);
delta=(2/3)*phi;

rbbw=size(bbw,1);
rbfr_s=size(bfr_s,1);
rbfr_w=size(bfr_w,1);
rbfr_wsid=size(bfr_wside,1);
L=2*rbbw+2*rbfr_s+2*rbfr_w+2*rbfr_wsid;

 pl=21;
 nd=size(ddata_s,1)+size(ddata_w,1);

tmtx=0:0.05:1;
% tmtx=1;
tm=size(tmtx,2);

for k=1:tm

 ah=0.1;
av=0.05;
he=H;
T=0.5;           % time period in second range 0.3-0.6 sec
coe=tmtx(k);
t=coe*T  ; 


%%%%%%%%%%%%%%%%%continuam equlibrium satisfied in triangles - soil %%%%%%%%%%%%%

aeq=sparse(2*E,3*N);
beq=sparse(2*E,1);

Ns=size(renamed_nodes_matrix_s,1);
Es=Ns/3;
aeq_s=sparse(2*Es,3*Ns);
beq_s=sparse(2*Es,1);
for i=1:Es

    %%%%%%%%%%% for soil %%%%%%%%%

%%% suction
y=centroids_s(i,2);
    al=0.1;
    n=5;
    q=0;
    ks=3e-5;
    gammaW=9.81;
    hw=0;
    Hw=hw*H;
    % Hw=-0.5;
    D=q/ks;
    Bo=exp(al*gammaW*(Hw-y));
    F=log(((1+D)*Bo)-D);
    if (Hw-y)>=0
        deno=1;
        Se=1;
    else
        deno=(1+(-F)^n).^((n-1)/n);
        Se=1/deno;
    end
    sig_s=(1/al)*(F/deno);
e_v=(G*gammaW/(gammaD))-1;    
gamma=(G+(e_v*Se))*gammaW/(1+e_v);


%%% overburden calcutn
z_vs=centroids_s(i,2);
ob=(gammaD*H)*((0.1074*(z_vs/H)^2)-(1.1421*z_vs/H)+1.0673);
p0=ob-sig_s;

% Vs varying %%
N60=((rad2deg(phi)-18.7)/(6.9*p0^(-0.12)))^(1/0.46);
Vs=85.88*N60^0.352;
Vp=1.87*Vs;

%%% dAMPING VARYING%%
Cu=2;
d50=0.5;
ro=0.5;
patm=101.325;

Di= (0.55*(Cu^0.1)*(d50^-0.3))*(-2.06*ro+2.43)*((p0/patm)^(0.72*ro-0.86));                % damping coeficient


ys1=2*pi*he/(T*Vs)*(((1+4*Di^2)^0.5+1)/(2+8*Di^2))^0.5;
ys2=-2*pi*he/(T*Vs)*(((1+4*Di^2)^0.5-1)/(2+8*Di^2))^0.5;


yp1=2*pi*he/(T*Vp)*(((1+4*Di^2)^0.5+1)/(2+8*Di^2))^0.5;
yp2=-2*pi*he/(T*Vp)*(((1+4*Di^2)^0.5-1)/(2+8*Di^2))^0.5;

cs=cos(ys1)*cosh(ys2);
ss=-sin(ys1)*sinh(ys2);

cp=cos(yp1)*cosh(yp2);
sp=-sin(yp1)*sinh(yp2);


    x1=renamed_nodes_matrix_s(T2_renamed(i,1),2);
    x2=renamed_nodes_matrix_s(T2_renamed(i,2),2);
    x3=renamed_nodes_matrix_s(T2_renamed(i,3),2);
    y1=renamed_nodes_matrix_s(T2_renamed(i,1),3);
    y2=renamed_nodes_matrix_s(T2_renamed(i,2),3);
    y3=renamed_nodes_matrix_s(T2_renamed(i,3),3);
    z_s=-centroids_s(i,2);

    eta1=y2-y3;     eta2=y3-y1;     eta3=y1-y2;
    zeta1=x3-x2;    zeta2=x1-x3;    zeta3=x2-x1;
    tw_area=(eta1*zeta2)-(eta2*zeta1);
    if tw_area<=0
        disp i
        disp tw_area
    end
   j=3*i-2;
   aeq_s((2*i-1):(2*i),(3*j-2):(3*j+6))=(1/tw_area)*[eta1 0   zeta1   eta2    0   zeta2   eta3    0   zeta3;
                                                    0 zeta1   eta1    0   zeta2   eta2    0   zeta3   eta3];


   

csz=cos(ys1*z_s/he)*cosh(ys2*z_s/he);
ssz=-sin(ys1*z_s/he)*sinh(ys2*z_s/he);
cpz=cos(yp1*z_s/he)*cosh(yp2*z_s/he);
spz=-sin(yp1*z_s/he)*sinh(yp2*z_s/he);
temp=ah/(cs^2+ss^2)*((cs*csz+ss*ssz)*cos(2*pi*t/T)+(ss*csz-cs*ssz)*sin(2*pi*t/T));
temp1=av/(cp^2+sp^2)*((cp*cpz+sp*spz)*cos(2*pi*t/T)+(sp*cpz-cp*spz)*sin(2*pi*t/T));

%  if(temp1<0)
%     temp1=-1*temp1;
% end   
    beq_s(2*i-1,1)=-temp*gamma;
    beq_s(2*i,1)=(1-temp1)*gamma;
     
end


%%%%%%%%%%%%%%%%%%%%% for wall %%%%%%%%%%%%%%%%%%%%%%%%%%

Nw=size(renamed_nodes_matrix_w,1);
Ew=Nw/3;
aeq_w=sparse(2*Ew,3*Nw);
beq_w=sparse(2*Ew,1);
for i=1:Ew

    %%%%%% for wall %%%%%%%%

Vs_w=2500;                % H/TVs=0.3  
Vp_w=3900;            % H/TVp=0.16


xi_w=0.05 ;                % damping coeficient


ys1_w=2*pi*he/(T*Vs_w)*(((1+4*xi_w^2)^0.5+1)/(2+8*xi_w^2))^0.5;
ys2_w=-2*pi*he/(T*Vs_w)*(((1+4*xi_w^2)^0.5-1)/(2+8*xi_w^2))^0.5;


yp1_w=2*pi*he/(T*Vp_w)*(((1+4*xi_w^2)^0.5+1)/(2+8*xi_w^2))^0.5;
yp2_w=-2*pi*he/(T*Vp_w)*(((1+4*xi_w^2)^0.5-1)/(2+8*xi_w^2))^0.5;

cs_w=cos(ys1_w)*cosh(ys2_w);
ss_w=-sin(ys1_w)*sinh(ys2_w);

cp_w=cos(yp1_w)*cosh(yp2_w);
sp_w=-sin(yp1_w)*sinh(yp2_w);


    x1=renamed_nodes_matrix_w(T1_renamed(i,1),2);
    x2=renamed_nodes_matrix_w(T1_renamed(i,2),2);
    x3=renamed_nodes_matrix_w(T1_renamed(i,3),2);
    y1=renamed_nodes_matrix_w(T1_renamed(i,1),3);
    y2=renamed_nodes_matrix_w(T1_renamed(i,2),3);
    y3=renamed_nodes_matrix_w(T1_renamed(i,3),3);
    z_w=-centroids_w(i,2);

    eta1=y2-y3;     eta2=y3-y1;     eta3=y1-y2;
    zeta1=x3-x2;    zeta2=x1-x3;    zeta3=x2-x1;
    tw_area=(eta1*zeta2)-(eta2*zeta1);
    if tw_area<=0
        disp i
        disp tw_area
    end
   j=3*i-2;
   aeq_w((2*i-1):(2*i),(3*j-2):(3*j+6))=(1/tw_area)*[eta1 0   zeta1   eta2    0   zeta2   eta3    0   zeta3;
                                                    0 zeta1   eta1    0   zeta2   eta2    0   zeta3   eta3];


csz_w=cos(ys1_w*z_w/he)*cosh(ys2_w*z_w/he);
ssz_w=-sin(ys1_w*z_w/he)*sinh(ys2_w*z_w/he);
cpz_w=cos(yp1_w*z_w/he)*cosh(yp2_w*z_w/he);
spz_w=-sin(yp1_w*z_w/he)*sinh(yp2_w*z_w/he);
temp=ah/(cs_w^2+ss_w^2)*((cs_w*csz_w+ss_w*ssz_w)*cos(2*pi*t/T)+(ss_w*csz_w-cs_w*ssz_w)*sin(2*pi*t/T));
temp1=av/(cp_w^2+sp_w^2)*((cp_w*cpz_w+sp_w*spz_w)*cos(2*pi*t/T)+(sp_w*cpz_w-cp_w*spz_w)*sin(2*pi*t/T));

%  if(temp1<0)
%     temp1=-1*temp1;
% end   
    beq_w(2*i-1,1)=-temp*gamma_w;
    beq_w(2*i,1)=(1-temp1)*gamma_w;
     
end
aeq(1:2*Es,1:3*Ns)=aeq_s;
aeq(2*Es+1:2*E,(3*Ns+1):3*N)=aeq_w;
beq(1:2*Es,1)=beq_s;
beq(2*Es+1:2*E,1)=beq_w;


%%%%%%%objective function%%%%%%%%%
f=zeros(1,6*N);
no_be=size(shared_edges,1);
for i=1:no_be

    lx_b=renamed_nodes_matrix_s(bbw(2*i-1,1),2)-renamed_nodes_matrix_s(bbw(2*i,1),2);
    ly_b=renamed_nodes_matrix_s(bbw(2*i-1,1),3)-renamed_nodes_matrix_s(bbw(2*i,1),3);
    % thb=atan(ly_b/lx_b);

% r=(lx_b^2+ly_b^2)^0.5;

    j=bbw(2*i-1,1);
    f(1,3*j-2)=0.5*ly_b;
if bbw(2*i-1,3)>=H2
f(1,3*j)=0.5*lx_b;
else
f(1,3*j)=-0.5*lx_b;
 end
    j=bbw(2*i,1);
    f(1,3*j-2)=0.5*ly_b;
if bbw(2*i,3)>=H2
f(1,3*j)=0.5*lx_b;
else
f(1,3*j)=-0.5*lx_b;
end
end

f=f.*1000;
accu=1e12;
      f=round(f.*accu)./accu;
      f=-f; 

      %%%%%%%% Discontinuity matrix  (soil)  %%%%%%%%%%
adis=sparse(4*nd,3*N);
bdis=sparse(4*nd,1);
nd_s=size(ddata_s,1);
adis_s=sparse(4*nd_s,3*Ns);
bdis_s=sparse(4*nd_s,1);
for i=1:nd_s
 
     thdis=ddata_s(i,5);
 
     T=[(sin(thdis))^2 , (cos(thdis))^2 , -sin(2*thdis); -0.5*sin(2*thdis) ,  0.5*sin(2*thdis), cos(2*thdis)];
 
            j=ddata_s(i,1);
 
             adis_s(4*i-3:4*i-2,3*j-2:3*j)=-T;
 
        j=ddata_s(i,2);
 
             adis_s(4*i-3:4*i-2,3*j-2:3*j)=T;
 
         j=ddata_s(i,3);
 
              adis_s(4*i-1:4*i,3*j-2:3*j)=-T;
 
        j=ddata_s(i,4);
 
             adis_s(4*i-1:4*i,3*j-2:3*j)=T;
end

%%%%%%% for wall %%%%%%%%%%%%%
nd_w=size(ddata_w,1);
adis_w=sparse(4*nd_w,3*Nw);
bdis_w=sparse(4*nd_w,1);
for i=1:nd_w
 
     thdis=ddata_w(i,5);
 
     T=[(sin(thdis))^2 , (cos(thdis))^2 , -sin(2*thdis); -0.5*sin(2*thdis) ,  0.5*sin(2*thdis), cos(2*thdis)];
 
            j=ddata_w(i,1);
 
             adis_w(4*i-3:4*i-2,3*j-2:3*j)=-T;
 
        j=ddata_w(i,2);
 
             adis_w(4*i-3:4*i-2,3*j-2:3*j)=T;
 
         j=ddata_w(i,3);
 
              adis_w(4*i-1:4*i,3*j-2:3*j)=-T;
 
        j=ddata_w(i,4);
 
             adis_w(4*i-1:4*i,3*j-2:3*j)=T;
end
adis(1:4*nd_s,1:3*Ns)=adis_s;
adis(4*nd_s+1:4*nd,(3*Ns+1):3*N)=adis_w;
bdis(1:4*nd_s,1)=bdis_s;
bdis(4*nd_s+1:4*nd,1)=bdis_w;


%%%%%%%%%%%%%   boundary conditions matrix %%%%%%%%%%%%%
abou=zeros(L,3*N);
bbou=zeros(L,1);

%%% at location of surcharge %%%%%

thb=beta;
qs=0;
 
  T2=[ -0.5*sin(2*thb) ,  0.5*sin(2*thb), cos(2*thb)];
  T1=[(sin(thb))^2 , (cos(thb))^2 , -sin(2*thb)];
i=1;
  for ib=1:size(bfr_s,1)
j=bfr_s(ib,1);
    abou(i,3*j-2:3*j)=T2;
    abou(i+1,3*j-2:3*j)=T1;
    bbou(i+1,1)=qs;
j=bfr_s(ib,2);
    abou(i+2,3*j-2:3*j)=T2;
    abou(i+3,3*j-2:3*j)=T1;
    bbou(i+3,1)=qs;

i=i+4;
  end

  %%%%% at top of wall %%%%

  for ib=1:size(bfr_w,1)
j=Ns+bfr_w(ib,1);
    abou(i,3*j)=1;
    abou(i+1,3*j-1)=1;
    bbou(i+1,1)=0;
j=Ns+bfr_w(ib,2);
    abou(i+2,3*j)=1;
    abou(i+3,3*j-1)=1;
    bbou(i+3,1)=0;

i=i+4;
  end

  %%%%%% at outer side of wall %%%%

  for ib=1:size(bfr_wside,1)
j=Ns+bfr_wside(ib,1);
    abou(i,3*j)=1;
    % abou(i+1,3*j-1)=1;
    abou(i+1,3*j-2)=1;
    bbou(i+1,1)=0;
j=Ns+bfr_wside(ib,2);
    abou(i+2,3*j)=1;
    % abou(i+3,3*j-1)=1;
    abou(i+3,3*j-2)=1;
    bbou(i+3,1)=0;

i=i+4;
  end

  %%% suction %%^%%%
  clear sig_s
for j=1:Ns
    y=renamed_nodes_matrix_s(j,3);
    al=0.1;
    n=5;
    q=0;
    ks=3e-5;
    gammaW=9.81;
    hw=0;
    Hw=hw*H;
     % Hw=-0.5;
    D=q/ks;
    Bo=exp(al*gammaW*(Hw-y));
    F=log(((1+D)*Bo)-D);
    if (Hw-y)>=0
        deno=1;
    else
        deno=(1+(-F)^n).^((n-1)/n);
    end
    sig_s(j,1)=(1/al)*(F/deno);
end

%%%%%% shear stress boundary conditions at the back of the wall %%%%%%%

  afin=zeros(4*size(shared_edges,1),3*N);
  bfin_l=zeros(4*size(shared_edges,1),1);
  bfin_u=zeros(4*size(shared_edges,1),1);

for  i=1:size(shared_edges,1)

   lx_b=renamed_nodes_matrix_s(bbw(2*i-1,1),2)-renamed_nodes_matrix_s(bbw(2*i,1),2);
   ly_b=renamed_nodes_matrix_s(bbw(2*i-1,1),3)-renamed_nodes_matrix_s(bbw(2*i,1),3);
   thb=atan(ly_b/lx_b);

T2=[-0.5*sin(2*thb) , 0.5*sin(2*thb), cos(2*thb)];
T1=[(sin(thb))^2 , (cos(thb))^2 , -sin(2*thb)];

j=bbw(2*i-1,1);

  afin(4*i-3,3*j-2:3*j)=T2+(T1)*tan(delta);
  bfin_l(4*i-3,1)=-inf;
  bfin_u(4*i-3,1)=(c*tan(delta)/tan(phi))-(sig_s(j,1)*tan(delta));
  afin(4*i-2,3*j-2:3*j)=-T2+(T1)*tan(delta);
  bfin_l(4*i-2,1)=-inf;
  bfin_u(4*i-2,1)=(c*tan(delta)/tan(phi))-(sig_s(j,1)*tan(delta));

    j=bbw(2*i,1);

  afin(4*i-1,3*j-2:3*j)=T2+(T1)*tan(delta);
  bfin_l(4*i-1,1)=-inf;
  bfin_u(4*i-1,1)=(c*tan(delta)/tan(phi))-(sig_s(j,1)*tan(delta));
  afin(4*i,3*j-2:3*j)=-T2+(T1)*tan(delta);
  bfin_l(4*i,1)=-inf;
  bfin_u(4*i,1)=(c*tan(delta)/tan(phi))-(sig_s(j,1)*tan(delta));

end

Aeq=[aeq;adis;abou;afin];
 
  Beq_L=[beq;bdis;bbou;bfin_l];
  Beq_U=[beq;bdis;bbou;bfin_u];


   Aeq=sparse(Aeq);
  Beq_L=sparse(Beq_L);
  Beq_U=sparse(Beq_U);

r=size(Aeq,1);

%%%%% Conic constraints %%%%%

asoc=sparse(3*N,3*N);
bsoc=sparse(3*N,1);


for i=1:N

    if i<=Ns
%%%%% for soil %%%%%
asoc(3*i-2,3*i-2)=sin(phi);
asoc(3*i-2,3*i-1)=sin(phi);
asoc(3*i-1,3*i-2)=-1;
asoc(3*i-1,3*i-1)=1;
asoc(3*i,3*i)=-2;
bsoc(3*i-2,1)=2*c*cos(phi)-(2*sig_s(i,1)*sin(phi));

elseif i>=Ns+1
     %%%%% for wall %%%%%
asoc(3*i-2,3*i-2)=sin(phi_w);
asoc(3*i-2,3*i-1)=sin(phi_w);
asoc(3*i-1,3*i-2)=-1;
asoc(3*i-1,3*i-1)=1;
asoc(3*i,3*i)=-2;
bsoc(3*i-2,1)=2*cw*cos(phi_w);
    end
end

Aeq=[Aeq;asoc];
Beq_L=[Beq_L;bsoc];
Beq_U=[Beq_U;bsoc];

  A2=[sparse(r,3*N);speye(3*N)];

  Aeq=[Aeq A2];

  Aeq=Aeq.*1000;
  Beq_L=Beq_L.*1000;
Beq_U=Beq_U.*1000;

accu=1e12;
Aeq=round(Aeq*accu)./accu;
Beq_L=round(Beq_L*accu)./accu;
Beq_U=round(Beq_U*accu)./accu;


%%%%%%%%%%%%%%5 socp optimization using mosek%%%%%%%%%%%%5

addpath 'C:\Program Files\Mosek\10.2\toolbox'\r2017aom\
clear prob;

[r, res] = mosekopt('symbcon');

prob.c   =  f;
prob.a   = Aeq;
prob.blc = Beq_L;
prob.buc = Beq_U;
prob.blx = -inf*ones(6*N,1);
prob.bux = inf*ones(6*N,1);

% Specify the cones as affine conic constraints.
temp=sparse(3*N,6*N);

pc=[res.symbcon.MSK_DOMAIN_QUADRATIC_CONE 3];

prob.accs=repmat(pc,1,N);

for i=1:N

temp(3*i-2,3*N+3*i-2)=1;
temp(3*i-1,3*N+3*i-1)=1;
temp(3*i,3*N+3*i)=1;

end

prob.f=sparse(temp);
[r,res]=mosekopt('minimize',prob);

% Display the primal solution.

x=res.sol.itr.xx;
PH(k)=res.sol.itr.pobjval/1000 ;

ka_h(k,1)=2*PH(k)/(gammaD*H^2);     % store Ka values at each t/T    




%%% resultant pressure calculation

P_h=zeros(size(bbw,1),1);
P_v=zeros(size(bbw,1),1);
ph_y=zeros(size(bbw,1),1);
Pa=zeros(size(bbw,1),3);
Pa_dist=zeros(size(bbw,1),3);
thta_b=zeros(no_be,1);

for i=1:no_be
    lx_b=abs(renamed_nodes_matrix_s(bbw(2*i-1,1),2)-renamed_nodes_matrix_s(bbw(2*i,1),2));
    ly_b=abs(renamed_nodes_matrix_s(bbw(2*i-1,1),3)-renamed_nodes_matrix_s(bbw(2*i,1),3));
r=(lx_b^2+ly_b^2)^0.5;
    j=bbw(2*i-1,1);
sigx_bnd=x(3*j-2,1);
sigy_bnd=x(3*j-1,1);
txy_bnd=x(3*j,1);
   if bbw(2*i-1,3)>=H2
       P_h(2*i-1,1)=(sigx_bnd*ly_b*0.5)+(txy_bnd*lx_b*0.5);
    P_v(2*i-1,1)=(sigy_bnd*lx_b*0.5)+(txy_bnd*ly_b*0.5);
    ph_y(2*i-1,1)=P_h(2*i-1,1)*renamed_nodes_matrix_s(bbw(2*i-1,1),3);
    ph_y1(2*i-1,1)=P_h(2*i-1,1)*renamed_nodes_matrix_s(bbw(2*i-1,1),3);
    Pa(2*i-1,1:2)=bbw(2*i-1,2:3);
    Pa(2*i-1,3)=P_h(2*i-1,1)/cos(theta1+delta);
    Pa_dist(2*i-1,1:2)=bbw(2*i-1,2:3);
    Pa_dist(2*i-1,3)=(Pa(2*i-1,3)/(0.5*r));
    % P_v(2*i-1,1)=P_h(2*i-1,1)*tan(theta1+delta);
   else
   
    P_h(2*i-1,1)=(sigx_bnd*ly_b*0.5)-(txy_bnd*lx_b*0.5);
    P_v(2*i-1,1)=-(sigy_bnd*lx_b*0.5)+(txy_bnd*ly_b*0.5);
    ph_y(2*i-1,1)=P_h(2*i-1,1)*renamed_nodes_matrix_s(bbw(2*i-1,1),3);
    ph_y2(2*i-1,1)=P_h(2*i-1,1)*renamed_nodes_matrix_s(bbw(2*i-1,1),3);
    Pa(2*i-1,1:2)=bbw(2*i-1,2:3);
    % Pa(2*i-1,3)=(P_h(2*i-1,1)/cos(theta2-delta));
    Pa(2*i-1,3)=(P_h(2*i-1,1)/cos(delta-theta2));
    Pa_dist(2*i-1,1:2)=bbw(2*i-1,2:3);
    Pa_dist(2*i-1,3)=(Pa(2*i-1,3)/(0.5*r));
    % P_v(2*i-1,1)=P_h(2*i-1,1)*tan(delta-theta2);
   end
   j=bbw(2*i,1);
sigx_bnd=x(3*j-2,1);
sigy_bnd=x(3*j-1,1);
txy_bnd=x(3*j,1);
   if bbw(2*i,3)>=H2
    
    P_h(2*i,1)=(sigx_bnd*ly_b*0.5)+(txy_bnd*lx_b*0.5);
    P_v(2*i,1)=(sigy_bnd*lx_b*0.5)+(txy_bnd*ly_b*0.5);
    ph_y(2*i,1)=P_h(2*i,1)*renamed_nodes_matrix_s(bbw(2*i,1),3);
    ph_y1(2*i,1)=P_h(2*i,1)*renamed_nodes_matrix_s(bbw(2*i,1),3);
    Pa(2*i,1:2)=bbw(2*i,2:3);
    Pa(2*i,3)=(P_h(2*i,1)/cos(theta1+delta));
    Pa_dist(2*i,1:2)=bbw(2*i,2:3);
    Pa_dist(2*i,3)=(Pa(2*i,3)/(0.5*r));
    % P_v(2*i,1)=P_h(2*i,1)*tan(theta1+delta);
   else

     P_h(2*i,1)=(sigx_bnd*ly_b*0.5)-(txy_bnd*lx_b*0.5);
     P_v(2*i,1)=-(sigy_bnd*lx_b*0.5)+(txy_bnd*ly_b*0.5);
     ph_y(2*i,1)=P_h(2*i,1)*renamed_nodes_matrix_s(bbw(2*i,1),3);
     ph_y2(2*i,1)=P_h(2*i,1)*renamed_nodes_matrix_s(bbw(2*i,1),3);
    Pa(2*i,1:2)=bbw(2*i,2:3);
    Pa(2*i,3)=(P_h(2*i,1)/cos(delta-theta2));
    Pa_dist(2*i,1:2)=bbw(2*i,2:3);
    Pa_dist(2*i,3)=(Pa(2*i,3)/(0.5*r));
    % P_v(2*i,1)=P_h(2*i,1)*tan(delta-theta2);
   end
end
   Ph_mtx=[bbw(:,2:3)  P_h];
Pv_mtx=[bbw(:,2:3)  P_v];
Pa_mtx=Pa;
Pa1h=[];  Pa1v=[];
Pa2h=[];  Pa2v=[];
Pa1=[];   Pa2=[];
for m=1:size(Ph_mtx,1)
    if Ph_mtx(m,2)>=H2
        Pa1h=[Pa1h;Ph_mtx(m,3)];
        Pa1v=[Pa1v;Pv_mtx(m,3)];
        Pa1=[Pa1;Pa_mtx(m,3)];
    else
        Pa2h=[Pa2h;Ph_mtx(m,3)];
        Pa2v=[Pa2v;Pv_mtx(m,3)];
        Pa2=[Pa2;Pa_mtx(m,3)];
    end
end

% store follwoing results at each t/T
P1(k,1)=sum(Pa1);
P2(k,1)=sum(Pa2);
Ph_t(k,1)=sum(Pa1h)+sum(Pa2h);
Pv_t(k,1)=sum(Pa1v)+sum(Pa2v);
Pr(k,1)=sqrt(Ph_t(k,1)^2+Pv_t(k,1)^2);
Kar(k,1)=2*Pr(k,1)/(gammaD*H^2);
theta_r2(k,1)=atan(Pv_t(k,1)/Ph_t(k,1));

%%%%% point of application

ybar=sum(ph_y)/sum(P_h);
zh(k,1)=ybar/H;
end

rslt(:,1)=ka_h  ; 
rslt(:,2)=P1;
rslt(:,3)=P2;
rslt(:,4)=Ph_t;
rslt(:,5)=Pv_t;
rslt(:,6)=Pr;
rslt(:,7)=Kar;
rslt(:,8)=theta_r2;
rslt(:,9)=zh;
rslt(:,10)=tmtx';

% max Ka out of the set corresponding to t/T=0 to 1 is identified and
% remaining result values are extracted correspondingly

    Ka_f=max(ka_h);
    [~,maxidx]=max(rslt(:,1));
    Pt=rslt(maxidx,2);
    Pb=rslt(maxidx,3);
    Ph=rslt(maxidx,4);
    Pv=rslt(maxidx,5);
    Pr=rslt(maxidx,6);
    Ka_r=rslt(maxidx,7);
    thta_r=rslt(maxidx,8);
    ZH=rslt(maxidx,9);
    tval=rslt(maxidx,10);