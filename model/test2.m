clc;clear all;close all

function qConfig = ikCodegen(endEffectorName,tform,weights,initialGuess)
    %#codegen
    persistent ik robot
    if isempty(robot)
        robot = loadrobot("kinovaGen3","DataFormat","row");
    end
    if isempty(ik)
        ik = inverseKinematics('RigidBodyTree',robot);
    end

    [qConfig,~] = ik(endEffectorName,tform,weights,initialGuess); 
end

robot = loadrobot("kinovaGen3",DataFormat="row");
showdetails(robot)

endEffectorName = "EndEffector_Link";
weights = [0.25 0.25 0.25 1 1 1];
initialGuess = [0 0 0 0 0 0 0];

targetPose = trvec2tform([0.35 -0.35 0]);
qConfig = ikCodegen(endEffectorName,targetPose,weights,initialGuess)

codegen ikCodegen -args {endEffectorName,targetPose,weights,initialGuess}

figure;
show(robot,qConfig);
hold on
plotTransforms(tform2trvec(targetPose),tform2quat(targetPose),FrameSize=0.5);
axis([-0.1 0.7 -0.5 0.5 -0.3 0.5])

%%

t = (0:0.2:10)'; % Time
count = length(t);
center = [0.3 0.3 0];
radius = 0.15;
theta = t*(2*pi/t(end));
points = center + radius*[cos(theta) sin(theta) zeros(size(theta))];

q0 = [0 0 0 0 0 0 0];
ndof = length(q0);
qs = zeros(count,ndof);
weights = [0 0 0 1 1 1];
endEffector = "EndEffector_Link";

qInitial = q0; % Use home configuration as the initial guess
tic
for i = 1:count
    % Solve for the configuration satisfying the desired end-effector
    % position
    point = points(i,:);
    qSol = ikCodegen(endEffector,trvec2tform(point),weights,qInitial);
    % Store the configuration
    qs(i,:) = qSol;
    % Start from prior solution
    qInitial = qSol;
end
loopTime = toc;
fprintf("Average (mexed) IK Execution Time: %f\n", loopTime/count);

robot = loadrobot("kinovaGen3",DataFormat="row");
% Show first solution and set view.
figure
show(robot,qs(1,:));
view(3)
ax = gca;
ax.Projection = "orthographic";
hold on
plot(points(:,1),points(:,2))
axis([-0.1 0.7 -0.5 0.5 -0.3 0.5])

% Iterate through the solutions.
framesPerSecond = 15;
r = rateControl(framesPerSecond);
for i = 1:count
    show(robot,qs(i,:),PreservePlot=false,FastUpdate=true);
    drawnow
    waitfor(r);
end