function R=ThVib_Group8_Rotation_function(NODE1,NODE2)
    X3=100;
    Y3=200;
    Z3=3000;
    %% Global axis
    X=[1 0 0];
    Y=[0 1 0];
    Z=[0 0 1];
    
    L=sqrt((NODE2(1)-NODE1(1))^2+(NODE2(2)-NODE1(2))^2+(NODE2(3)-NODE1(3))^2);

    coord_global=[NODE1(1) NODE1(2) NODE1(3); NODE2(1) NODE2(2) NODE2(3);X3 Y3 Z3]; %Node 1 Node 2 and the third point--> element 1
    %x,y,z local axis
    x=[(NODE2(1)-NODE1(1))/L,(NODE2(2)-NODE1(2))/L,(NODE2(3)-NODE1(3))/L];
    d2=[NODE2(1)-NODE1(1),NODE2(2)-NODE1(2),NODE2(3)-NODE1(3)];
    d3=[X3-NODE1(1),Y3-NODE1(2),Z3-NODE1(3)];
    y=cross(d2,d3)/norm(cross(d2,d3));
    z=cross(x,y);
    R=[dot(X,x) dot(Y,x) dot(Z,x);dot(X,y) dot(Y,y) dot(Z,y);dot(X,z) dot(Y,z) dot(Z,z)];
end