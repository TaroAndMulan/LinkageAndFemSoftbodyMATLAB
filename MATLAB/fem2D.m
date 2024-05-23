function linkages(scene)
if nargin < 1
    scene = 0;
end

pins = [];
particles = [];
switch scene
    case 0
        % Crank-rocker
        % Bottom link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; % display vertices
        % Left link
        links(2).angle = pi/2;
        links(2).pos = [-1;0];
        links(2).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]';
        % Right link
        links(3).angle = pi/2;
        links(3).pos = [1;0];
        links(3).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]';
        % Top link
        links(4).angle = 0;
        links(4).pos = [-1;1];
        links(4).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        
        % Which link is grounded?
        grounded = 1;
        % Which link is the driver?
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % Bottom-left
        pins(1).links = [1,2];
        pins(1).pts = [0,0;0,0]';
        % Bottom-right
        pins(2).links = [1,3];
        pins(2).pts = [2,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [1,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [1+rand(1),0;2,0]'; % pin location on link3 is randomized
        
        % List of tracer particles for display
        particles(1).link = 4; % which link?
        particles(1).pt = [0.5;0.1]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        particles(2).link = 4;
        particles(2).pt = [2.5;-0.1];
        particles(2).ptsWorld = zeros(2,0);
    case 1

        % Drag-link

        % s
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [0;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]'; % display vertices
        % q
        links(2).angle = pi/2;
        links(2).pos = [0;0];
        links(2).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        % l
        links(3).angle = pi/2;
        links(3).pos = [1;0];
        links(3).verts = [0,-0.1;4,-0.1;4,0.1;0,0.1]';
        % p
        links(4).angle = 0;
        links(4).pos = [0;3];
        links(4).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        
        % Which link is grounded?
        grounded = 1;
        % Which link is the driver?
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % Bottom-left
        pins(1).links = [1,2];
        pins(1).pts = [0,0;0,0]';
        % Bottom-right
        pins(2).links = [1,3];
        pins(2).pts = [1,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [3,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [4,0;3,0]'; % pin location on link3 is randomized
        
        % List of tracer particles for display
        particles(1).link = 4; % which link?
        particles(1).pt = [0;0.1]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        particles(2).link = 4;
        particles(2).pt = [3;-0.1];
        particles(2).ptsWorld = zeros(2,0);
    case 2
        % Double-rocker

        % Bottom link q
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [0;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;2,-0.1;0,0.1;2,0.1]'; % display vertices
        % Left link l
        links(2).angle = pi/3;
        links(2).pos = [0;0];
        links(2).verts = [0,-0.1;5,-0.1;5,0.1;0,0.1]';
        % Right link p 
        links(3).angle = pi/2;
        links(3).pos = [2;0];
        links(3).verts = [0,-0.1;3.5,-0.1;3.5,0.1;0,0.1]';
        % Top link s
        links(4).angle = 0;
        links(4).pos = [2.5;4.33];
        links(4).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]';
        
        % Which link is grounded?
        grounded = 1;
        % Which link is the driver?
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % Bottom-left
        pins(1).links = [1,2];
        pins(1).pts = [0,0;0,0]';
        % Bottom-right
        pins(2).links = [2,4];
        pins(2).pts = [5,0;0,0]';
        % Left-top
        pins(3).links = [1,3];
        pins(3).pts = [2,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [3.5,0;1,0]'; % pin location on link3 is randomized
        
        % List of tracer particles for display
        particles(1).link = 4; % which link?
        particles(1).pt = [0;0]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        particles(2).link = 4;
        particles(2).pt = [1;0];
        particles(2).ptsWorld = zeros(2,0);
    case 3
        % Hoekens

        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [0;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;2,-0.1;0,0.1;2,0.1]'; % display vertices

        links(2).angle = pi/2;
        links(2).pos = [0;0];
        links(2).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]';
   

        links(3).angle = pi/2;
        links(3).pos = [2;0];
        links(3).verts = [0,-0.1;2.5,-0.1;2.5,0.1;0,0.1]';

        links(4).angle = 0;
        links(4).pos = [0;1];
        links(4).verts = [0,-0.1;5,-0.1;5,0.1;0,0.1]';
        

        grounded = 1;
 
        driver = 2;
        
        pins(1).links = [1,2];
        pins(1).pts = [0,0;0,0]';

        pins(2).links = [1,3];
        pins(2).pts = [2,0;0,0]';

        pins(3).links = [2,4];
        pins(3).pts = [1,0;0,0]';

        pins(4).links = [3,4];
        pins(4).pts = [2.5,0;2.5,0]';
        
        % List of tracer particles for display
        particles(1).link = 4; 
        particles(1).pt = [5;0.1]; 
        particles(1).ptsWorld = zeros(2,0); 

    case 4
        % Peaucellier-Lipkin

        % grounded
        links(1).angle = 0;
        links(1).pos = [0;0]; 
        links(1).verts = [0,0]'; 
        
        links(2).angle = 0;
        links(2).pos = [1;0]; 
        links(2).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]'; 
        
        links(3).angle = pi/4;
        links(3).pos = [2;0]; 
        links(3).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; 

        links(4).angle = -pi/4;
        links(4).pos = [2;0]; 
        links(4).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; 

        links(5).angle = pi/8;
        links(5).pos = [2;0]; 
        links(5).verts = [0,-0.1;3.7,-0.1;3.7,0.1;0,0.1]';  

        links(6).angle = -pi/8;
        links(6).pos = [2;0]; 
        links(6).verts = [0,-0.1;3.7,-0.1;3.7,0.1;0,0.1]'; 
        
        links(7).angle = -pi/4;
        links(7).pos = [3.414;1.414]; 
        links(7).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; 

        links(8).angle = pi/4;
        links(8).pos = [3.414;-1.414]; 
        links(8).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; 

        grounded = 1;

        driver = 2;
   
        pins(1).links = [1,2];
        pins(1).pts = [1,0;0,0]';

        pins(2).links = [2,3];
        pins(2).pts = [1,0;0,0]';

        pins(3).links = [2,4];
        pins(3).pts = [1,0;0,0]';

        pins(4).links = [1,5];
        pins(4).pts = [0,0;0,0]';

        pins(5).links = [1,6];
        pins(5).pts = [0,0;0,0]';

        pins(6).links = [3,5];
        pins(6).pts = [2,0;3.7,0]';

        pins(7).links = [4,6];
        pins(7).pts = [2,0;3.7,0]';

        pins(8).links = [3,7];
        pins(8).pts = [2,0;0,0]';

        pins(9).links = [4,8];
        pins(9).pts = [2,0;0,0]';

        pins(10).links = [7,8];
        pins(10).pts = [2,0;2,0]';


        % List of tracer particles for display
        particles(1).link = 7; 
        particles(1).pt = [2;0];
        particles(1).ptsWorld = zeros(2,0); 

    case 5
        % Klann

        links(1).angle = 0.29;
        links(1).pos = [3.4;8.03]; 
        links(1).verts = [0,-0.1;6.888,-0.1;6.36857,7.356232]';

        links(2).angle = 2*pi/3;
        links(2).pos = [10;10]; 
        links(2).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';

        links(3).angle = 2*pi/3;
        links(3).pos = [3.4;8.03]; 
        links(3).verts = [0,-0.1;3.59,-0.1;3.59,0.1;0,0.1]';

        links(4).angle = pi;
        links(4).pos = [5;11]; 
        links(4).verts = [0,-0.1;6.56,0.53;12.3636,-0.1;12.3636,0.1;6.56,0.73;0,0.1]';

        links(5).angle = 2*pi/3;
        links(5).pos = [7.4;16.9]; 
        links(5).verts = [0,-0.1;5.79,-0.1;5.79,0.1;0,0.1]';

        links(6).angle = 3*pi/2;
        links(6).pos = [3;18]; 
        links(6).verts = [0,-0.1;10.04,-0.1;10.04,0.1;0,0.1]';

        links(7).angle = 1.4*pi;
        links(7).pos = [-5;11]; 
        links(7).verts = [0,-0.1;10.04,-0.1;10.04,0.1;0,0.1]';

        grounded = 1;

        driver = 2;

        pins(1).links = [1,2];
        pins(1).pts = [6.8877,0;0,0]';

        pins(7).links = [1,5];
        pins(7).pts = [6.36857,7.356232;0,0]';

        pins(4).links = [1,3];
        pins(4).pts = [0,0;0,0]';
        
        pins(2).links = [2,4];
        pins(2).pts = [3,0;0,0]';

        pins(3).links = [3,4];
        pins(3).pts = [3.59,0;6.56,0.73]';

        pins(6).links = [5,6];
        pins(6).pts = [5.79,0;0,0]';

        pins(5).links = [4,6];
        pins(5).pts = [12.3636,0;10.04,0]';

        pins(8).links = [6,7];
        pins(8).pts = [10.04,0;0,0]';

        particles(1).link = 7; 
        particles(1).pt = [10.04;0.1]; 
        particles(1).ptsWorld = zeros(2,0); 
 

    case 10
end

% Initialize
for i = 1 : length(links)
    links(i).grounded = (i == grounded);
    links(i).driver = (i == driver);
    % These target quantities are only used for grounded and driver links
    links(i).angleTarget = links(i).angle;
    links(i).posTarget = links(i).pos;
end
drawScene(links,pins,particles,true);

% lsqnonlin options
if verLessThan('matlab','8.1')
    opt = optimset(...
        'Jacobian','on',...
        'DerivativeCheck','off',...
        'Display','off'); % final-detailed iter-detailed off
else
    opt = optimoptions('lsqnonlin',...
        'Jacobian','on',...
        'DerivativeCheck','off',...
        'Display','off'); % final-detailed iter-detailed off
end

% Simulation loop
t = 0; % current time
T = 1; % final time
dt = 0.01; % time step
angVel = 2*pi; % driver angular velocity

while t < T
    %links(driver).angleTarget = links(driver).angleTarget + dt*angVel;
    switch scene
        case 0
            links(driver).angleTarget = links(driver).angleTarget + dt*angVel;
        case 1
            links(driver).angleTarget = links(driver).angleTarget + dt*angVel;
        case 2
            min = pi/8;
            max = pi/3;
            f=1;
            links(driver).angleTarget = 0.5*(max-min)*sin(2*pi*f*t+1.5*pi)+0.5*(max+min);
        case 3
            links(driver).angleTarget = links(driver).angleTarget + 2*dt*angVel;
        case 4
            min = -pi/3;
            max = pi/3;
            f=2;
            links(driver).angleTarget = 0.5*(max-min)*sin(2*pi*f*t)+0.5*(max+min);
        case 5
            links(driver).angleTarget = links(driver).angleTarget - 4*dt*angVel;
    end

    % Solve for linkage orientations and positions
    [links,feasible] = solveLinkage(links,pins,opt);
    % Update particle positions
    particles = updateParticles(links,particles);
    % Draw scene
    drawScene(links,pins,particles);
    % Quit if over-constrained
    if ~feasible
        break;
    end
    t = t + dt;
end

end

%%
function [R,dR] = rotationMatrix(angle)
c = cos(angle);
s = sin(angle);
% Rotation matrix
R = zeros(2);
R(1,1) = c;
R(1,2) = -s;
R(2,1) = s;
R(2,2) = c;
if nargout >= 2
    % Rotation matrix derivative
    dR = zeros(2);
    dR(1,1) = -s;
    dR(1,2) = -c;
    dR(2,1) = c;
    dR(2,2) = -s;
end
end

%%
function [links,feasible] = solveLinkage(links,pins,opt)
nlinks = length(links);
% Extract the current angles and positions into a vector
angPos0 = zeros(3*nlinks,1);
for i = 1 : nlinks
    link = links(i);
    ii = (i-1)*3+(1:3);
    angPos0(ii(1)) = link.angle;
    angPos0(ii(2:3)) = link.pos;
end
% Limits
lb = -inf(size(angPos0));
ub =  inf(size(angPos0));
% Solve for angles and positions
[angPos,r2] = lsqnonlin(@(angPos)objFun(angPos,links,pins),angPos0,lb,ub,opt);
% If the mechanism is feasible, then the residual should be zero
feasible = true;
if r2 > 1e-6
    fprintf('Mechanism is over constrained!\n');
    feasible = false;
end
% Extract the angles and positions from the values in the vector
for i = 1 : length(links)
    ii = (i-1)*3+(1:3);
    links(i).angle = angPos(ii(1));
    links(i).pos = angPos(ii(2:3));
end
end

%%
function [f,J] = objFun(angPos,links,pins)
nlinks = length(links);
npins = length(pins);
% Temporarily change angles and positions of the links. These changes will
% be undone when exiting this function.
for i = 1 : nlinks
    ii = (i-1)*3+(1:3);
    links(i).angle = angPos(ii(1));
    links(i).pos = angPos(ii(2:3));
end

% Evaluate constraints
ndof = 3*nlinks;
ncon = 3 + 3 + 2*npins; % 3 for ground, 3 for driver, 2*npins for pins
f = zeros(ncon,1);
J = zeros(ncon,ndof);
k = 0;
% Some angles and positions are fixed
for i = 1 : nlinks
    link = links(i);
    if link.grounded || link.driver
        % Grounded and driver links have their angles and positions
        % prescribed.
        f(k+1,    1) = link.angle - link.angleTarget;
        f(k+(2:3),1) = link.pos - link.posTarget;
        % The Jacobian of this constraint is the identity matrix
        J(k+(1:3),k+(1:3)) = eye(3);
        k = k + 3;
    end
end
% Pin constraints
for i = 1 : npins
    pin = pins(i);
    rows = k+(1:2); % row index of this pin constraint
    indLinkA = pin.links(1); % array index of link A
    indLinkB = pin.links(2); % array index of link B
    linkA = links(indLinkA);
    linkB = links(indLinkB);
    [Ra,dRa] = rotationMatrix(linkA.angle);
    [Rb,dRb] = rotationMatrix(linkB.angle);
    % Local positions
    ra = pin.pts(:,1);
    rb = pin.pts(:,2);
    % World positions
    xa = Ra * ra + linkA.pos;
    xb = Rb * rb + linkB.pos;
    p = xa(1:2) - xb(1:2);
    f(rows,1) = p;
    % Column indices for the angles and positions of links A and B
    colAngA = (indLinkA-1)*3 + 1;
    colPosA = (indLinkA-1)*3 + (2:3);
    colAngB = (indLinkB-1)*3 + 1;
    colPosB = (indLinkB-1)*3 + (2:3);
    % The Jacobian of this constraint is the partial derivative of f wrt
    % the angles and positions of the two links.
    J(rows,colAngA) = dRa * ra;
    J(rows,colPosA) = eye(2);
    J(rows,colAngB) = -dRb * rb;
    J(rows,colPosB) = -eye(2);
    k = k + 2;
end
end

%%
function particles = updateParticles(links,particles)
% Transform particle position from local to world
for i = 1 : length(particles)
    particle = particles(i);
    link = links(particle.link);
    R = rotationMatrix(link.angle);
    x = R * particle.pt + link.pos;
    % Append world position to the array (grows indefinitely)
    particles(i).ptsWorld(:,end+1) = x;
end
end

%%
function drawScene(links,pins,particles,initialize)
if nargin < 4
    initialize = false;
end
if initialize
    clf;
else
    cla;
end
hold on;
grid on;
% Draw links
for i = 1 : length(links)
    link = links(i);
    R = rotationMatrix(link.angle);
    % Draw frame
    p = link.pos; % frame origin
    s = 0.1; % frame display size
    px = p + s*R(:,1); % frame x-axis
    py = p + s*R(:,2); % frame y-axis
    plot([p(1),px(1)],[p(2),px(2)],'r','LineWidth',3);
    plot([p(1),py(1)],[p(2),py(2)],'g','LineWidth',3);
    % Draw link geometry
    if link.grounded
        color = [1 0 0];
    elseif link.driver
        color = [0 1 0];
    else
        color = [0 0 1];
    end
    E = [R,link.pos;0,0,1]; % transformation matrix
    vertsLocal = [link.verts;ones(1,size(link.verts,2))];
    vertsWorld = E * vertsLocal;
    plot(vertsWorld(1,[1:end,1]),vertsWorld(2,[1:end,1]),'Color',color);
end
% Draw pins
for i = 1 : length(pins)
    pin = pins(i);
    linkA = links(pin.links(1));
    linkB = links(pin.links(2));
    Ra = rotationMatrix(linkA.angle);
    Rb = rotationMatrix(linkB.angle);
    xa = Ra * pin.pts(:,1) + linkA.pos;
    xb = Rb * pin.pts(:,2) + linkB.pos;
    plot(xa(1),xa(2),'co','MarkerSize',10,'MarkerFaceColor','c');
    plot(xb(1),xb(2),'mx','MarkerSize',10,'LineWidth',2);
end
% Draw particles
for i = 1 : length(particles)
    particle = particles(i);
    if ~isempty(particle.ptsWorld)
        plot(particle.ptsWorld(1,:),particle.ptsWorld(2,:),'k');
        plot(particle.ptsWorld(1,end),particle.ptsWorld(2,end),'ko');
    end
end
axis equal;
drawnow;
end
