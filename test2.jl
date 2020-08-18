using LibFEM,Plots
# Units system mm/tonne/s/K
# ===============================================================================
# 2D FRAME
# ===============================================================================
# (2)---(3)
#      /
#     /
#    /
#   /
#  /
# /
#(1)
# ===============================================================================
# The element used here is a linear 2D straight truss with constant cross section
# The degrees of freedom are the u & v displacements
# ===============================================================================
#PARAMETERS
L=10000.; #length in mm
A=100.; #cross-sectional area in mm^2
E=210000.; #modulus of elasticity in MPa [steel]
FM=10000.; #force modulus in N
# ===============================================PRE-PROCESSING==================
# DEFINING AND DISCRETIZING[MESHING] THE STRUCTURE
# connectivity table
# elt||node_i||node_j
# 1|1|3
# 2|2|3
# grid()
X1pos=0.;Y1pos=0.
X2pos=0.;Y2pos=L
X3pos=L;Y3pos=L
#lengths
L1=D2_TrussElementLength(X1pos,Y1pos,X3pos,Y3pos); #length of element 1
L2=D2_TrussElementLength(X2pos,Y2pos,X3pos,Y3pos); #length of element 2
#APPLYING GEOMETRIC&MATERIAL PROPERTIES
A1=sqrt(2)*A; #cross-sectional area of element 1
A2=A; #cross-sectional area of element 2
E1=E; #material of element 1
E2=E; #material of element 2
#writing-defining the element stiffness matrices
K1=D2_TrussElementStiffness(E1,A1,L1,45);
println("\nK1=\r")
display(K1)
K2=D2_TrussElementStiffness(E2,A2,L2,0);
println("\nK2=\r")
display(K2)
#ASSEMBLING THE GLOBAL STIFFNESS MATRIX
#matrices initialization
K=zeros(6,6);K1P=zeros(6,6);
#positionning stiffness matrices
K1P=D2_TrussAssemble(K,K1,1,3)
println("\nK1P=\r")
display(K1P)
K=zeros(6,6);K2P=zeros(6,6);
K2P=D2_TrussAssemble(K,K2,2,3)
println("\nK2P=\r")
display(K2P)
#assembling
K=K1P+K2P
println("\nK=\r")
display(K)
#SOLVING DISPACEMENT EQUATIONS
#extracting displacement submatrix via index vector
K_s=K[5:6,5:6]
#Setting-up the force subvector by applying Load & Boundary Conditions[LBC]]
F_s=[0 FM]'
#solving by gaussian elimination
U_s=K_s\F_s
#SOLVING FORCE EQUATIONS
#setting-up the global nodal displacement vector
U=[0 0 0 0 U_s[1:2]']'
println("\nU=\r")
display(U)
#computing the global nodal force vector
F=K*U
println("\nF=\r")
display(F)
#COMPUTING STRESSES
#writing the element nodal displacement vectors
U1=[U[1] U[2] U[5] U[6]]'
U2=[U[3] U[4] U[5] U[6]]'
#computing element strains
ϵ1=D2_TrussElementStrain(L1,45,U1)
println("\nϵ1=\r")
display(ϵ1)
ϵ2=D2_TrussElementStrain(L2,0,U2)
println("\nϵ2=\r")
display(ϵ2)#computing element forces
f1=D2_TrussElementForce(E1,A1,L1,45,U1)
f2=D2_TrussElementForce(E2,A2,L2,0,U2)
#computing element stresses
σ1=D2_TrussElementStress(E1,L1,45,U1)
println("\nσ1=\r")
display(σ1)
σ2=D2_TrussElementStress(E2,L2,0,U2)
println("\nσ2=\r")
display(σ2)
