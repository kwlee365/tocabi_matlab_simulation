tocabi = importrobot('dyros_tocabi.urdf');

tocabi.DataFormat = 'column';  % 또는 'row'로 설정 가능 (q 형식과 일치시킬 것)
tocabi.Gravity = [0; 0; -9.81];  % 중력 설정

q = zeros(40,1);      % 관절 위치

q(4) = 1.0;

q(7) = 0.7;

M = massMatrix(tocabi, q)
gtau = gravityTorque(tocabi,q)
eeName = 'L_Foot_Link'; 
J = geometricJacobian(tocabi, q, eeName)


show(tocabi, q)
axis equal

showdetails(tocabi)