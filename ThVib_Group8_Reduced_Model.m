clear all
close all

%% PARAMETERS TO CHOOSE:

Nelement=4; %Nb of elements for each beams

NB_EIGENVALUES=8; %Nb of eigenvalue taking into consideration

NB_MODES=10;%Only for Craig & Bampton

%% 
%% I-Modelling and modal analysis
%% 
%% 1/PARAMETERS TOWER


rho=8*10^3;        % [kg/m3] material density
E=210e9;% [Pa] Young's modulus
nu=0.3;%Poisson coefficient
Jx= 2.25*(0.07/2)^4;%Moment of torsion
Iy=0.07*0.07^3/12; %Moment of inertia
Iz=0.07*0.07^3/12;%Moment of inertia
A=0.07*0.07; % [m2] Cross section area
rt=sqrt(Jx/A);%Gyroscopic radius
G=E/(2*(1+nu)); %IN Pa shear modulus of steel

%% 2/PARAMETERS FE
% * List of coordinates of all the beams: [NODES1 0 NODES2]

COORD_BEAMS=10^-3*[0 0 0 0 750 750 5000;10000 0 0 0 9250 750 5000;10000 10000 0 0 9250 9250 5000;0 10000 0 0 750 9250 5000;750 750 5000 0 9250 750 5000;9250 750 5000 0 9250 9250 5000;9250 9250 5000 0 750 9250 5000;750 9250 5000 0 750 750 5000;0 0 0 0 9250 750 5000;10000 0 0 0 9250 9250 5000;10000 10000 0 0 750 9250 5000;0 10000 0 0 750 750 5000;1500 1500 10000 0 8500 1500 10000;8500 1500 10000 0 8500 8500 10000;8500 8500 10000 0 1500 8500 10000;1500 8500 10000 0 1500 1500 10000;750 750 5000 0 1500 1500 10000;9250 750 5000 0 8500 1500 10000;9250 9250 5000 0 8500 8500 10000;750 9250 5000 0 1500 8500 10000;9250 750 5000 0 1500 1500 10000;9250 9250 5000 0 8500 1500 10000;750 9250 5000 0 8500 8500 10000;750 750 5000 0 1500 8500 10000;2250 2250 15000 0 7750 2250 15000;7750 2250 15000 0 7750 7750 15000;7750 7750 15000 0 2250 7750 15000;2250 7750 15000 0 2250 2250 15000;1500 1500 10000 0 2250 2250 15000;8500 1500 10000 0 7750 2250 15000;8500 8500 10000 0 7750 7750 15000;1500 8500 10000 0 2250 7750 15000;1500 1500 10000 0 7750 2250 15000;8500 1500 10000 0 7750 7750 15000;8500 8500 10000 0 2250 7750 15000;1500 8500 10000 0 2250 2250 15000;3000 3000 20000 0 7000 3000 20000;7000 3000 20000 0 7000 7000 20000;7000 7000 20000 0 3000 7000 20000;3000 7000 20000 0 3000 3000 20000;2250 2250 15000 0 3000 3000 20000;7750 2250 15000 0 7000 3000 20000;7750 7750 15000 0 7000 7000 20000;2250 7750 15000 0 3000 7000 20000;7750 2250 15000 0 3000 3000 20000;7750 7750 15000 0 7000 3000 20000;2250 7750 15000 0 7000 7000 20000;2250 2250 15000 0 3000 7000 20000;3000 3000 20000 0 3000 3000 23000;7000 3000 20000 0 7000 3000 23000;7000 7000 20000 0 7000 7000 23000;3000 7000 20000 0 3000 7000 23000;3000 3000 23000 0 7000 3000 23000;7000 3000 23000 0 7000 7000 23000;7000 7000 23000 0 3000 7000 23000;3000 7000 23000 0 3000 3000 23000;3000 3000 20000 0 5000 -3000 20000;7000 3000 20000 0 5000 -3000 20000;3000 3000 23000 0 5000 -3000 20000;7000 3000 23000 0 5000 -3000 20000;7000 7000 20000 0 5000 13000 20000;3000 7000 20000 0 5000 13000 20000;7000 7000 23000 0 5000 13000 20000;3000 7000 23000 0 5000 13000 20000];

%% 
% * List of coordinates of all NODES: [NODES]

NODES=[COORD_BEAMS(:,1:3);COORD_BEAMS(:,5:7)];
NODES=unique(NODES,'rows'); %To keep only one time one node (because one node can appear several times because it can belongs to 1,2 or 3 beams)
nb_nodes_tower=length(NODES);
nb_beams_tower=length(COORD_BEAMS);

%% Reduction
%%

index_condensed=1;
index_remaining=1;

for i = 1:length(NODES)
    if NODES(i,3)<20-10^-4;
        condensed_nodes(index_condensed, :)=NODES(i,:);
        index_condensed=index_condensed+1;
    else
          remaining_nodes(index_remaining, :) = NODES(i,:);
          index_remaining=index_remaining+1;
     end
end 
NODES(:,1:3)=[remaining_nodes;condensed_nodes];


%% 3/Discretisation of the tower
%%
% * Parameters of the discretization: We have to test differents values of Nelement 
% in order to discuss on this parameter
% 
% What I think:
% 
% -If Nelement is to small (2,3,4) , the discrezation is too large and so 
% the results will not be correct (far from the reality)
% 
% -If Nelement too big: it will require a long time to compute the matix 
% and it is too much --> problem of convergence
% 
% 
% 
% To check in running NX for differents values of discretization.
% 
% * Discritization through x-axis in the local axis
% 
% |-----------------------|-----------------------|-----------------------|-----------------------|
% 
% N1                      N1,1                    N1,2                   N1,3                     N2             
%                       
% 
% N1,1, N1,2, N1,3 : secoundary nodes

for i=1:nb_beams_tower
    R=ThVib_Group8_Rotation_function(COORD_BEAMS(i,1:3),COORD_BEAMS(i,5:7));%Rotation matrix to go from global to local
    COORD_NODE1=R*COORD_BEAMS(i,1:3)';%Coordinate of node 1 in local axis
    COORD_NODE2=R*COORD_BEAMS(i,5:7)'; %Coordinate of node 1 in local axis
    X_DISCRETISATION=linspace(COORD_NODE1(1), COORD_NODE2(1), Nelement); %Discretization throug x-axis
    
    counter_secondary_nodes=2;%We start in taking the second node because the first one is a N1 a primary nodes
    counter_colomn=1;
    while counter_secondary_nodes<Nelement
        NODES_FEM_BIS(i,counter_colomn:counter_colomn+2)=inv(R)*[X_DISCRETISATION(counter_secondary_nodes),COORD_NODE1(2),COORD_NODE1(3)]';
        %[X_DISCRETISATION(counter_colomn),COORD_NODE1(2),COORD_NODE1(3)]=[x-coordinate given by linspace, y-coordinate: same than the primary nodes, z-coordinate: same than the primary nodes]
        %y and z are constant in the local axis.
        %*inv(R) --> to go back into the global axis
        counter_colomn=counter_colomn+3;%Increase by three to fill the next three coloms
        counter_secondary_nodes=counter_secondary_nodes+1;%Increase by one to take the next x-coordinate of the next secondary node
    end 

end

%% 
% Above, in NODES_FEM_BIS, all the secondary nodes are refered: in one row, 
% all the nodes belonging to a same beam are refered.
% 
% 
% 
% -Get the list of coordinates of each nodes after applying FEM
% 
% Below, the purpose of the code is to reference all the node in the list 
% called NODES_FEM: one row=one node

counter_beam=1;%we start in considering the first beam
counter_row=1;%We start to fill the first row
while counter_beam<=length(NODES_FEM_BIS) %We want to go through all the beam
    for k=1:3:(Nelement-2)*3 %k is the index of the colomn
        NODES_FEM(counter_row,1:3)=NODES_FEM_BIS(counter_beam,k:k+2);%We reference each secondary nodes in one big matrix
        counter_row=counter_row+1;%One you fill the row with one node, we increase by one to fill the next row
    end
    counter_beam=counter_beam+1;%Once we reference all the secondary nodes of the beam we go on the next one

end

NODES_FEM=[NODES;NODES_FEM];%We add to the secondary nodes, the primary nodes
%NODES_FEM=[PRIMARY NODES;SECONDARY NODES]

%% Get the list of coordinate for each elements after applying FEM


NODES_FEM_TER=[COORD_BEAMS(:,1:3),NODES_FEM_BIS,COORD_BEAMS(:,5:7)];
%NODES_FEM_TER=[Node 1, secondary nodes, Node 2]
counter_row=1;%We start to fill the first colomn
counter_beam=1;% We start taking the first beam
while counter_beam<=length(NODES_FEM_TER) %We want to go through all the beams (64)

    for k=1:3:Nelement*3-3 %k represents the index of the colomn
        %We know that two consecutives nodes in the matrix are linked and so they form one element
        COORD_FEM(counter_row,1:3)=NODES_FEM_TER(counter_beam,k:k+2);%Nodes 1
        
        COORD_FEM(counter_row,4:6)=NODES_FEM_TER(counter_beam,k+3:k+5);%Nodes 2
        
         counter_row=counter_row+1;%One you have [Nodes 1 Nodes2], we increase the counter_row to fill the next row with the next element
    end
    counter_beam=counter_beam+1;%One you went through all the elements of one beam, you can go through the next beam
end

        
%%  Coordinates for the reduction method

condensed_nodes=[condensed_nodes;NODES_FEM(length(NODES)+1:end,:)];
%% Plot after applying FEM

figure(1)
%Plot the segment betwenn two nodes
for i=1:length(COORD_BEAMS)
    hold on
 
    plot3([COORD_BEAMS(i,1),COORD_BEAMS(i,5)],[COORD_BEAMS(i,2),COORD_BEAMS(i,6)],[COORD_BEAMS(i,3),COORD_BEAMS(i,7)],'r')
    
end

%Plot the nodes
for i=1:length(NODES_FEM)
    hold on
    
    scatter3(NODES_FEM(i,1),NODES_FEM(i,2),NODES_FEM(i,3),'o')
    
end
     

%% Localization matrix
%  Degrees of freedom

NB_NODES_FEM=length(NODES_FEM);
NB_DOF=6*NB_NODES_FEM;%6 DOF for each nodes

% List of degrees of freedom (DoFs related to each node)
dofList = [(1:6:NB_DOF-5)' (2:6:NB_DOF-4)' (3:6:NB_DOF-3)' (4:6:NB_DOF-2)' (5:6:NB_DOF-1)' (6:6:NB_DOF)'];

% * Localisation matrix: [degrees of freedom of node 1, degrees of freedom of 
% node 2]= 12 colomns=6*2nodes
% 
% For each element, we are going to try to find the degrees of freedom of 
% both nodes and then we can put together those DOF using locel(i,:)= [dofs1,dofs2].
% 
% To do that, for each element, we are going to:
% 
%     -Take the coordinate of node 1 and node2
% 
%     -Then we went to determine the degrees of freedom of those nodes. So 
% we are going to try to indentify those nodes in the list of node: NODES_FEM.
% 
%     -One you have the row index of Node 1 and Node 2, you can have the 
% DOF of those nodes going in dofList(index1,:) and dofList(index2,:).

for i=1:length(COORD_FEM)
    NODE1=COORD_FEM(i,1:3);
    NODE2=COORD_FEM(i,4:6);
    index1=1;%We start at the first row of NODES_FEM
    while (NODES_FEM(index1,1)~=NODE1(1)) || (NODES_FEM(index1,2)~=NODE1(2)) || (NODES_FEM(index1,3)~=NODE1(3))
       %While you do not have identify your node in the list of nodes, you will increanse the index, to try the next values
       % | : or --> While x or y or z are not equal we increase the index to go further
        index1=index1+1;      
    
    end
    %Exactly the same but for the seconde node
    index2=1;
    while (NODES_FEM(index2,1)~=NODE2(1)) ||(NODES_FEM(index2,2)~=NODE2(2))||(NODES_FEM(index2,3)~=NODE2(3))
        index2=index2+1;      
    end
    INDEX_COORD(i,:)=[index1 index2];
    dofs1 = dofList(index1, :); %Degrees of freedom of the first node
    dofs2 = dofList(index2, :); %Degrees of freedom of the seconde node
    locel(i,:) = [dofs1,dofs2];
end

    
%% Stiffeness and mass matrix

 Ks=zeros(NB_DOF, NB_DOF);%Initial matrix
 Ms=zeros(NB_DOF, NB_DOF);%Initial matrix

for i=1:length(COORD_FEM) %Going through all the elements
    
    L=sqrt((COORD_FEM(i,1)-COORD_FEM(i,4))^2+(COORD_FEM(i,2)-COORD_FEM(i,5))^2+(COORD_FEM(i,3)-COORD_FEM(i,6))^2);
    l=L; %Length of your element
    
    %Elementary stiffeness matrix in LOCAL AXIS !!
    Kel =[E*A/l 0 0 0 0 0 -E*A/l 0 0 0 0 0;
    0 12*E*Iz/l^3 0 0 0 6*E*Iz/l^2 0 -12*E*Iz/l^3 0 0 0 6*E*Iz/l^2;
    0 0 12*E*Iy/l^3 0 -6*E*Iy/l^2 0 0 0 -12*E*Iy/l^3 0 -6*E*Iy/l^2 0;
    0 0 0 G*Jx/l 0 0 0 0 0 -G*Jx/l 0 0;
    0 0 -6*E*Iy/l^2 0 4*E*Iy/l 0 0 0 6*E*Iy/l^2 0 2*E*Iy/l 0;
    0 6*E*Iz/l^2 0 0 0 4*E*Iz/l 0 -6*E*Iz/l^2 0 0 0 2*E*Iz/l;
    -E*A/l 0 0 0 0 0 E*A/l 0 0 0 0 0;
    0 -12*E*Iz/l^3 0 0 0 -6*E*Iz/l^2 0 12*E*Iz/l^3 0 0 0 -6*E*Iz/l^2;
    0 0 -12*E*Iy/l^3 0 6*E*Iy/l^2 0 0 0 12*E*Iy/l^3 0 6*E*Iy/l^2 0;
    0 0 0 -G*Jx/l 0 0 0 0 0 G*Jx/l 0 0;
    0 0 -6*E*Iy/l^2 0 2*E*Iy/l 0 0 0 6*E*Iy/l^2 0 4*E*Iy/l 0;
    0 6*E*Iz/l^2 0 0 0 2*E*Iz/l 0 -6*E*Iz/l^2 0 0 0 4*E*Iz/l];
    

 %Elementary mass matrix in LOCAL AXIS !!
    Mel = rho*A*l*...
    [1/3 0 0 0 0 0 1/6 0 0 0 0 0;
    0 13/35 0 0 0 11*l/210 0 9/70 0 0 0 -13*l/420;
    0 0 13/35 0 -11*l/210 0 0 0 9/70 0 13*l/420 0;
    0 0 0 rt^2/3 0 0 0 0 0 rt^2/6 0 0;
    0 0 -11*l/210 0 l^2/105 0 0 0 -13*l/420 0 -l^2/140 0;
    0 11*l/210 0 0 0 l^2/105 0 13*l/420 0 0 0 -l^2/140;
    1/6 0 0 0 0 0 1/3 0 0 0 0 0;
    0 9/70 0 0 0 13*l/420 0 13/35 0 0 0 -11*l/210;
    0 0 9/70 0 -13*l/420 0 0 0 13/35 0 11*l/210 0;
    0 0 0 rt^2/6 0 0 0 0 0 rt^2/3 0 0;
    0 0 13*l/420 0 -l^2/140 0 0 0 11*l/210 0 l^2/105 0;
    0 -13*l/420 0 0 0 -l^2/140 0 -11*l/210 0 0 0 l^2/105];

    % ROtation matrix
    R=ThVib_Group8_Rotation_function(COORD_FEM(i,1:3),COORD_FEM(i,4:6));
    
    T=blkdiag(R,R,R,R);%Block diagonal matrix --> 12*12 matrix
    %Here, we need a 12*12 rotation matrix to do the matricial product with Kel or Mel
    %One block R is for 1 block (x,y,z) either for the translation or the rotation
%% Stiffness and mass matrices of the element expressed in structural axes
    KeS=T'*Kel*T;
    MeS=T'*Mel*T;
    
    %Here we use the locel matrix to insert the elementary matrix at the good place in the big matrix
    Ks(locel(i,:),locel(i,:)) = Ks(locel(i,:),locel(i,:)) + KeS;
    
    Ms(locel(i,:),locel(i,:)) = Ms(locel(i,:),locel(i,:)) + MeS;
 
end

%% Find the nodes at both extremity
indexExtremity=zeros(2,1);
i=1;
index=1;
while index<=2
    if NODES_FEM(i,2)==-3
        indexExtremity(index,:)=i;
        index=index+1;
        
    elseif NODES_FEM(i,2)==13
        indexExtremity(index,:)=i;
        index=index+1;
    end
    i=i+1;
    
end
        
%% Lumped mass
%ADD the lumped mass at both extremity of the tower: 13,14 are the index of the nodes at both extremities
Ms(dofList(indexExtremity(1),1:3),dofList(indexExtremity(1),1:3))=Ms(dofList(indexExtremity(1),1:3),dofList(indexExtremity(1),1:3))+300;
%Ks(dofList(13,1:3),dofList(13,1:3))=Ks(dofList(13,1:3),dofList(13,1:3))+300;
Ms(dofList(indexExtremity(2),1:3),dofList(indexExtremity(2),1:3))=Ms(dofList(indexExtremity(2),1:3),dofList(indexExtremity(2),1:3))+300;
%Ks(dofList(14,1:3),dofList(14,1:3))=Ks(dofList(14,1:3),dofList(14,1:3))+300;



%% Boundary conditions
% Below, the purpose of the code is to find the nodes wich are clamped in order 
% to delete those elements in the main matrix Ks and Ms

i=1;
index=1;%We are going through all the nodes to find the clamped nodes
NODES_CLAMPED=zeros(4,1);
while index<5
    if NODES_FEM(i,3)==0 %One node is clamped if the z-coordinate is null
        NODES_CLAMPED(index)=i;%We add the index of each clamped nodes
        index=index+1;
    end
    i=i+1;
end

% %For each clamped nodes, we have to remove the rows and colomns related to the degrees of freedom o those nodes

for i=4:-1:1

    DOF_CLAMPED=dofList(NODES_CLAMPED,:);

    Ks(DOF_CLAMPED(i,6:-1:1),:)=[];
    Ks(:,DOF_CLAMPED(i,6:-1:1))=[];

    Ms(DOF_CLAMPED(i,6:-1:1),:)=[];
    Ms(:,DOF_CLAMPED(i,6:-1:1))=[];
end



   
NDOF=length(Ks);%Number of degrees of freedom after removing the clamped nodes
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%...............REDUCTION...................%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REDUCTION GUYAN-IRON
%%

KRR_GI=Ks(1:6*length(remaining_nodes),1:6*length(remaining_nodes));%Remaining part
KCC_GI=Ks(6*length(remaining_nodes)+1:end,6*length(remaining_nodes)+1:end);%Condensed part
KRC_GI=Ks(1:6*length(remaining_nodes),6*length(remaining_nodes)+1:end);
KCR_GI=KRC_GI';


RCR_GI=-inv(KCC_GI)*KCR_GI;
R_GI=[eye(length(Ks)-length(RCR_GI));RCR_GI];%Reduction matrix for Guyan-Iron

Ks_reduced_GI=R_GI'*Ks*R_GI;%Ks reduction --> decrease the number of rows/columns in Ks 
Ms_reduced_GI=R_GI'*Ms*R_GI;%Ms reduction --> decrease the number of rows/columns in Ms


%% REDUCTION CRAIG & BAMPTON
%%

KRR_CB=Ks(1:6*length(remaining_nodes),1:6*length(remaining_nodes));%Remaining part
KCC_CB=Ks(6*length(remaining_nodes)+1:end,6*length(remaining_nodes)+1:end);%condenced part

MRR_CB=Ms(1:6*length(remaining_nodes),1:6*length(remaining_nodes));%Remaining part
MCC_CB=Ms(6*length(remaining_nodes)+1:end,6*length(remaining_nodes)+1:end);%Condensed part

KRC_CB=Ks(1:6*length(remaining_nodes),6*length(remaining_nodes)+1:end);
KCR_CB=KRC_CB';

MRC_CB=Ms(1:6*length(remaining_nodes),6*length(remaining_nodes)+1:end);
MCR_CB=MRC_CB';

RCR_CB=-inv(KCC_CB)*KCR_CB;


%% SOLVE THE PROBLEM

%% Solve the exact solution of the problem: without reduction

[MODE, D]=eigs(Ks, Ms,NB_EIGENVALUES,'sm'); %Solve the problem
%eigs(Ks, Ms,8,'sm'): 8: ask for the 8 first eigenvalues and 'sm': to get a rank list

w = sqrt(diag(D));  % [rad/s]
f = w /2/pi;        % [Hz]

%% 

% Solution for the condensed nodes: in order to compute the phi, we have 
% to compute some modes of the static part (condensed part) and then add it to 
% the reduction matrix

[MODE_condensed_CB, D_condensed_CB]=eigs(KCC_CB, MCC_CB,NB_MODES,'sm');


%% 
% Solution for Guyan Iron

[MODE, D]=eigs(Ks, Ms,NB_EIGENVALUES,'sm'); %Solve the problem
[MODE_reduction_GI, D_reduction_GI]=eigs(Ks_reduced_GI, Ms_reduced_GI,NB_EIGENVALUES,'sm');%Solving again the problem but with the reduced matrix
%eigs(Ks, Ms,8,'sm'): 8: ask for the 8 first eigenvalues and 'sm': to get a rank list

w_reduction_GI = sqrt(diag(D_reduction_GI));  %
f_reduction_GI = w_reduction_GI /2/pi;        %Frequencies after reduction

eigenvalues_errors_GI=(f_reduction_GI-f)./f;

%% 
% Solution for Graig&Bampton

size_eye=size(RCR_CB);
R_CRAIG_BAMPTON=[eye(size_eye(2)),zeros(size_eye(2),NB_MODES);RCR_CB,MODE_condensed_CB];%Reduction matrix for CB (same than GI but with the eigenvectors of the condensed part in addition)

%Reduced matrix
Ks_reduced_CB=R_CRAIG_BAMPTON'*[KRR_CB,KRC_CB;KCR_CB,KCC_CB]*R_CRAIG_BAMPTON;
Ms_reduced_CB=R_CRAIG_BAMPTON'*[MRR_CB,MRC_CB;MCR_CB,MCC_CB]*R_CRAIG_BAMPTON;

%Solve the problem with the reduced matrix
[MODE_reduction_CB, D_reduction_CB]=eigs(Ks_reduced_CB, Ms_reduced_CB,NB_EIGENVALUES,'sm');
w_reduced_CB = sqrt(diag(D_reduction_CB));  % [rad/s]
f_reduced_CB = w_reduced_CB /2/pi; %reduced frequencies for CB


eigenvalues_errors_CB=(f_reduced_CB-f)./f;

%% 
% Damping

DAMPING_COEFFICIENT=[0.01,0.01]; %Damping for the two first modes
W_CB=2*pi*f_reduced_CB; %Eight first pulsations for CB (undamped problem)
W_GI=2*pi*f_reduction_GI;%Eight first pulsations for GI (undamped problem)
%% 
% * Calculate the coefficents a and b

a_CB=2*DAMPING_COEFFICIENT(1)/(W_CB(1)+W_CB(2));
b_CB=a_CB*(W_CB(1)*W_CB(2));

a_GI=2*DAMPING_COEFFICIENT(1)/(W_GI(1)+W_GI(2));
b_GI=a_GI*(W_GI(1)*W_GI(2));

C_reduced_GI=a_GI*Ks_reduced_GI+b_GI*Ms_reduced_GI;
C_reduced_CB=a_CB*Ks_reduced_CB+b_CB*Ms_reduced_CB;
%% 
% * Other damping ratio

for i=3:length(W_CB)
    DAMPING_COEFFICIENT_CB(i)=0.5*(a_CB*W_CB(i)+b_CB/W_CB(i));
end
for i=3:length(W_GI)
    DAMPING_COEFFICIENT_GI(i)=0.5*(a_GI*W_GI(i)+b_GI/W_GI(i));
end
%%

%% 
% * Load force

temporal_step=0.01; %Time step
T=30;%Final time
t=0:temporal_step:T;
  
  period=0.4;
force=0.4 * 10^3; %Force 0.4 kN
%Initialization of the load vectors
P_GI=zeros(length(Ks_reduced_GI),length(t));
P_CB=zeros(length(Ks_reduced_CB),length(t));

%Harmonic excitation for CB
P_CB(dofList(indexExtremity(1),3),:)=force*sin(2*pi/period*t);%-2 because the index changed after having removed the DOF of the clamped nodes
P_CB(dofList(indexExtremity(2),3),:)=-force*sin(2*pi/period*t);

%Harmonic excitation for GI
P_GI(dofList(indexExtremity(1),3),:)=force*sin(2*pi/period*t);%-2 because the index changed after having removed the DOF of the clamped nodes
P_GI(dofList(indexExtremity(2),3),:)=-force*sin(2*pi/period*t);


%% DIRECT TIME INTEGRATION
%%

% Guyan-Iron
        %Newmark reduction

gamma=0.5;
beta=0.25;

q0_GI=zeros(length(Ks_reduced_GI),1); %Initialization off the displacement at 0
dq0_GI=zeros(length(Ks_reduced_GI),1);%Initialization off the velocity at 0


d2q0_GI=inv(Ms_reduced_GI)*(P_GI(:,1)-C_reduced_GI*dq0_GI-Ks_reduced_GI*q0_GI);%Initialization of q dot dot 0

d2q_GI(1,:)=d2q0_GI;
dq_GI(1,:)=dq0_GI;
q_GI(1,:)=q0_GI;

d2q_star_GI(1,:)=d2q0_GI;
dq_star_GI(1,:)=dq0_GI;
q_star_GI(1,:)=q0_GI;

S_GI=Ms_reduced_GI+gamma*temporal_step*C_reduced_GI+beta*temporal_step^2*Ks_reduced_GI;
S_GI=inv(S_GI);
for temp=1:length(t)-1
    
    %Prediction
    dq_star_GI(temp+1,:)=dq_GI(temp,:)+(1-gamma)*temporal_step*d2q_GI(temp,:);
    q_star_GI(temp+1,:)=q_GI(temp,:)+temporal_step*dq_GI(temp,:)+(0.5-beta)*temporal_step^2*d2q_GI(temp,:);
    %COmputation of accelerations
    d2q_GI(temp+1,:)=S_GI*(P_GI(:,temp+1)-C_reduced_GI*dq_star_GI(temp+1,:)'-Ks_reduced_GI*q_star_GI(temp+1,:)');
    %Correction
    dq_GI(temp+1,:)=dq_star_GI(temp+1,:)+temporal_step*gamma*d2q_GI(temp+1,:);
    q_GI(temp+1,:)=q_star_GI(temp+1,:)+temporal_step^2*beta*d2q_GI(temp+1,:);
end
    
%% 
% Craig and Bampton

        %Newmark reduction

gamma=0.5;
beta=0.25;

q0_CB=zeros(length(Ks_reduced_CB),1); %Initialization off the displacement at 0
dq0_CB=zeros(length(Ks_reduced_CB),1);%Initialization off the velocity at 0


d2q0_CB=inv(Ms_reduced_CB)*(P_CB(:,1)-C_reduced_CB*dq0_CB-Ks_reduced_CB*q0_CB);%Initialization of q dot dot 0

d2q_CB(1,:)=d2q0_CB;
dq_CB(1,:)=dq0_CB;
q_CB(1,:)=q0_CB;

d2q_star_CB(1,:)=d2q0_CB;
dq_star_CB(1,:)=dq0_CB;
q_star_CB(1,:)=q0_CB;

S_CB=Ms_reduced_CB+gamma*temporal_step*C_reduced_CB+beta*temporal_step^2*Ks_reduced_CB;
S_CB=inv(S_CB);
for temp=1:length(t)-1
    
    %Prediction
    dq_star_CB(temp+1,:)=dq_CB(temp,:)+(1-gamma)*temporal_step*d2q_CB(temp,:);
    q_star_CB(temp+1,:)=q_CB(temp,:)+temporal_step*dq_CB(temp,:)+(0.5-beta)*temporal_step^2*d2q_CB(temp,:);
    %COmputation of accelerations
    d2q_CB(temp+1,:)=S_CB*(P_CB(:,temp+1)-C_reduced_CB*dq_star_CB(temp+1,:)'-Ks_reduced_CB*q_star_CB(temp+1,:)');
    %Correction
    dq_CB(temp+1,:)=dq_star_CB(temp+1,:)+temporal_step*gamma*d2q_CB(temp+1,:);
    q_CB(temp+1,:)=q_star_CB(temp+1,:)+temporal_step^2*beta*d2q_CB(temp+1,:);
end
%% 
% Guyan-Iron

 %%Dynamic response -- X axis
figure

plot(t,q_GI(:,dofList(indexExtremity(1),1)),'g')
legend('Time integration')
title('Guyan-Iron: Dynamic response X-axis')

%%Dynamic response -- Y axis
figure
plot(t,q_GI(:,dofList(indexExtremity(1),2)),'g')
legend('Time integration')
title('Guyan-Iron: Dynamic response Y-axis')

%%Dynamic response -- Z axis
figure
plot(t,q_GI(:,dofList(indexExtremity(1),3)),'g')
legend('Time integration')
title('Guyan-Iron: Dynamic response Z-axis')
  
%% 
% CRAIG&BAMPTON

   %%Dynamic response -- X axis
figure

plot(t,q_CB(:,dofList(indexExtremity(1),1)),'g')
legend('Time integration')
title('CRAIG&BAMPTON: Dynamic response X-axis')
%%Dynamic response -- Y axis
figure

plot(t,q_CB(:,dofList(indexExtremity(1),2)),'g')
legend('Time integration')
title('CRAIG&BAMPTON:Dynamic response Y-axis')
%%Dynamic response -- Z axis
figure

plot(t,q_CB(:,dofList(indexExtremity(1),3)),'g')
legend('Time integration')
title('CRAIG&BAMPTON: Dynamic response Z-axis')
  


