function [ R ] = Q2R( q )
q0=q(1);qx=q(2);qy=q(3);qz=q(4);
R=[q0^2+qx^2-qy^2-qz^2,2*(qx*qy-q0*qz),2*(qx*qz+q0*qy);...
    2*(qy*qx+q0*qz),q0^2-qx^2+qy^2-qz^2,2*(qy*qz-q0*qx);...
    2*(qz*qx-q0*qy),2*(qz*qy+q0*qx),q0^2-qx^2-qy^2+qz^2];
end

