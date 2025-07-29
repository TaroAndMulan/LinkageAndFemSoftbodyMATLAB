function fem2D(arg)

%create square mesh
if (arg==0 || arg==1)
[tri, v] = createSquareMesh(10,10, 1,1);
end
if (arg==2)
[tri, v] = createCircleMesh();
end

trimesh(tri, v(:,1), v(:,2));
numVerts = size(v,1);
pos = zeros(numVerts,2);

%floor is at zero
odeFun = @(time, state)(femOde(time,state,tri,v,pos, arg));
outputFun = @(time, state, flag)(femOutputFcn(time,state,flag,pos, tri, v));
options = odeset('OutputFcn', outputFun);

y0 = [reshape(v', 2*numVerts,1); zeros(2*numVerts,1)];
[time, state] = ode45(odeFun, [0, 10], y0, options);

pos(1:numVerts,:) = [state(end, 1:2:2*numVerts)' state(end, 2:2:2*numVerts)'];
trimesh(tri, pos(:,1), pos(:,2));
    
end

function y = femOde(time, state, tri, v, pos, arg)
    numVerts = size(v,1);
    pos(1:numVerts,:) = [state(1:2:2*numVerts) state(2:2:2*numVerts)];
    
    %hack to fix object to the floor
    %just set y positions of vertices to zero
    floorVert = abs(v(:,2)) < 1e-6; %some arbitrary toleranc
    %floorVert = abs(pos(:,2))< 1;
    %set positions to zero
    pos(floorVert,2) = 0;
    
    %set velocities to zero
    vel = [state((2*numVerts+1):2:end) state((2*numVerts+2):2:end)];
    vel(floorVert,2) = 0;
    y = [reshape(vel', 2*numVerts,1) ; femAccelerations(tri, v, pos, arg)];
    
end

%compute FEM accelerations
function a = femAccelerations(tri, v, pos, arg)

numTris = size(tri, 1);
numVerts = size(v,1);
f = zeros(numVerts,2);
m = zeros(numVerts, 2);
density = 1;
g = [0 -9.8]';

for i=1:numTris
    
    %Edge vectors
    E1 = (v(tri(i,2), :)-v(tri(i,1), :))';
    E2 = (v(tri(i,3), :)-v(tri(i,1), :))';
    E3 = (v(tri(i,3), :)-v(tri(i,2), :))';
    e1 = (pos(tri(i,2), :)-pos(tri(i,1), :))';
    e2 = (pos(tri(i,3), :)-pos(tri(i,1), :))';
    e3 = (pos(tri(i,3), :)-pos(tri(i,2), :))';
    area = 0.5*det([e1 e2]);
    me = density*area;
    
    %compute F= eE^-1
    F = [e1 e2]/([E1 E2]);
   
    %edit the cauchyStress method to add material models
    
    if (arg==0)
        stress = cauchyStress(F);
    end
    if (arg==1)
        stress = VenantKirchoffStress(F);
    end

    %compute forces for each edge
    
    %edge 1
    n1 = [e1(2) ; -e1(1)];
    
    %Set fe1 to the force on the first edge of the triangle (from node 1 to
    %node 2)
    fe1 = (stress*norm(e1))*(n1./norm(n1)); 
    
    
    f(tri(i,2), :) = f(tri(i,2), :) - 0.5*fe1';
    f(tri(i,1), :) = f(tri(i,1), :) - 0.5*fe1';
    
    %edge2
    n2 = [e3(2) ; -e3(1)];
    
    %Set fe2 to the force on the first edge of the triangle (from node 2 to
    %node 3)
    fe2 = (stress*norm(e3))*(n2./norm(n2)); 
    f(tri(i,3), :) = f(tri(i,3), :) - 0.5*fe2';
    f(tri(i,2), :) = f(tri(i,2), :) - 0.5*fe2';
    
    
    n3 = [-e2(2) ; e2(1)];
    
    %Set fe3 to the force on the first edge of the triangle (from node 3 to
    %node 1)
    fe3 = (stress*norm(e2))*(n3./norm(n3)); 
    f(tri(i,3), :) = f(tri(i,3), :) - 0.5*fe3';
    f(tri(i,1), :) = f(tri(i,1), :) - 0.5*fe3';
    
    
    %distribute mass to all vertices
    m(tri(i,1), : ) = m(tri(i,1), : ) + [me me]/3;
    m(tri(i,2), : ) = m(tri(i,2), : ) + [me me]/3;
    m(tri(i,3), : ) = m(tri(i,3), : ) + [me me]/3;
   
    
end

%reshape force vector and compute accelerations, add gravity here
f = reshape(f', 2*numVerts,1);
m = reshape(m', 2*numVerts,1);
a = f./m + repmat(g, numVerts,1);

end

function status = femOutputFcn(time, state, flag, pos, tri, v)

if strcmp(flag, 'done') == 0
numVerts = size(v,1);
hold on
clf
pos(1:numVerts,:) = [state(1:2:2*numVerts, end) state(2:2:2*numVerts, end)];
trimesh(tri, pos(:,1), pos(:,2));

hold off
xlim([-10, 10]);
ylim([0, 10]);
drawnow
end

status = 0;
end



function stress = cauchyStress(F)
    
    %Set P = to the Piola-Kirchoff 1 stress for a neohookean model.
    %Use parameters mu = 50, lambda = 50 to test
    mu = 50;
    lambda = 50;
    J = det(F);
    P = mu*(F-inv(transpose(F)))+(lambda*log(J)*inv(transpose(F)));
    stress = (P*F')./J;
    
end

function stress = VenantKirchoffStress(F)
    
    %Set P = to the Piola-Kirchoff 1 stress for a neohookean model.
    %Use parameters mu = 50, lambda = 50 to test
    E = 0.5*(transpose(F)*F-eye(2));
    mu = 50;
    lambda = 50;
    J = det(F);
    P = F*(2*mu*E+lambda*trace(E)*eye(2));
    stress = (P*F')./J;
    
end



function [tri, v] = createSquareMesh(width, height, dx, dy)

[X,Y] = meshgrid(1:dx:width, 1:dy:height);

v = [X(:) Y(:)];
tri = delaunay(v(:,1), v(:,2));
numVerts = size(v,1);
%center mesh
v(:,1) = v(:,1) - sum(v(:,1))./numVerts;
v(:,2) = v(:,2) - min(v(:,2));
%v(:,2) = v(:,2) + 10*ones(numVerts,1);

end

function [tri, v] = createCircleMesh()


M = 4 ;
N = 10 ;
R = 1
nR = linspace(0,R,M) ;
nT = linspace(0,2*pi,N) ;
%nT = pi/180*(0:NT:theta) ;
[R, T] = meshgrid(nR,nT) ;
% Convert grid to cartesian coordintes
X = R.*cos(T); 
Y = R.*sin(T);

v = [X(:) Y(:)];

tri = delaunay(v(:,1), v(:,2));
numVerts = size(v,1);
%center mesh
v(:,1) = v(:,1) - sum(v(:,1))./numVerts;
v(:,2) = v(:,2) - min(v(:,2));
end

