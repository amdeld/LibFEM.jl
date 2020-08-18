module LibFEM
function D1_SpringAssemble(K,k,i,j)
   #function prototype: D1_SpringAssemble(K,k,i,j)
   #This function assembles the element stiffness
   #matrix k of the spring with nodes i & j into the
   #global stiffness matrix K.
   #This function returns the global stiffness matrix K 
   #after the element stiffness matrix k is assembled.
   K[i,i] = K[i,i] + k[1,1]
   K[i,j] = K[i,j] + k[1,2]
   K[j,i] = K[j,i] + k[2,1]
   K[j,j] = K[j,j] + k[2,2]
   return K   
   end
   export D1_SpringAssemble
   function D1_SpringElementForce(k,u)
      #function prototype: D1_SpringElementForce(k,u)
      #This function returns the element nodal force
      #vector given the element stiffness matrix k 
      #and the element nodal displacement vector u.
      return k * u
      end
      export D1_SpringElementForce
function D2_TrussAssemble(K,k,i,j)
   #function prototype : D2_TrussAssemble(K,k,i,j)
   #This function assembles the element stiffness
   #matrix k of the plane truss element with nodes
   #i & j into the global stiffness matrix K.
   #This function returns the global stiffness
   #matrix K after the element stiffness matrix
   #k is assembled.
   K[2*i-1,2*i-1] = K[2*i-1,2*i-1] + k[1,1]
   K[2*i-1,2*i] = K[2*i-1,2*i] + k[1,2]
   K[2*i-1,2*j-1] = K[2*i-1,2*j-1] + k[1,3]
   K[2*i-1,2*j] = K[2*i-1,2*j] + k[1,4]
   K[2*i,2*i-1] = K[2*i,2*i-1] + k[2,1]
   K[2*i,2*i] = K[2*i,2*i] + k[2,2]
   K[2*i,2*j-1] = K[2*i,2*j-1] + k[2,3]
   K[2*i,2*j] = K[2*i,2*j] + k[2,4]
   K[2*j-1,2*i-1] = K[2*j-1,2*i-1] + k[3,1]
   K[2*j-1,2*i] = K[2*j-1,2*i] + k[3,2]
   K[2*j-1,2*j-1] = K[2*j-1,2*j-1] + k[3,3]
   K[2*j-1,2*j] = K[2*j-1,2*j] + k[3,4]
   K[2*j,2*i-1] = K[2*j,2*i-1] + k[4,1]
   K[2*j,2*i] = K[2*j,2*i] + k[4,2]
   K[2*j,2*j-1] = K[2*j,2*j-1] + k[4,3]
   K[2*j,2*j] = K[2*j,2*j] + k[4,4]
   return K
   end
   export D2_TrussAssemble
   function D2_TrussElementForce(E,A,L,theta,u)
#function prototype: D2_TrussElementForce(E,A,L,theta,u)
#This function returns the element force
#given the modulus of elasticity E; the
#cross-sectional area A; the length L;
#the angle theta [in degrees], & the
#element nodal displacement vector u.
x = theta * pi/180
C = cos(x)
S = sin(x)
return E*A/L*[-C -S C S]* u
end
export D2_TrussElementForce
function D2_TrussElementLength(x1,y1,x2,y2)
#function prototype: D2_TrussElementLength(x1,y1,x2,y2)
#This function returns the length of the
#plane truss element whose first node has
#coordinates [x1,y1] & second node has
#coordinates [x2,y2].
return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1))
end
export D2_TrussElementLength
function D2_TrussElementStiffness(E,A,L,theta)
#function prototype: D2_TrussElementStiffness(E,A,L,theta)
#This function returns the element
#stiffness matrix for a plane truss
#element with modulus of elasticity E;
#cross-sectional area A; length L; &
#angle theta [in degrees].
#The size of the element stiffness
#matrix is 4 x 4.
x = theta*pi/180
C = cos(x)
S = sin(x)
return E*A/L*[C*C C*S -C*C -C*S ; C*S S*S -C*S -S*S
   -C*C -C*S C*C C*S ; -C*S -S*S C*S S*S]
end
export D2_TrussElementStiffness
function D2_TrussElementStrain(L,theta,u)
#function prototype D2_TrussElementStrain(L,theta,u)
#This function returns the element strain
#given the length L,the angle theta [in degrees]
#and the element nodal displacement vector u.
x = theta * pi/180
C = cos(x)
S = sin(x)
return 1/L*[-C -S C S]* u
end
export D2_TrussElementStrain
function D2_TrussElementStress(E,L,theta,u)
#function prototype: D2_TrussElementStress(E,L,theta,u)
#This function returns the element stress
#given the modulus of elasticity E; the
#the length L; the angle theta (in
#degrees); & the element nodal
#displacement vector u.
x = theta * pi/180
C = cos(x)
S = sin(x)
return E/L*[-C -S C S]* u
end
export D2_TrussElementStress
end # module
