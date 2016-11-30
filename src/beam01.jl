function beam01(E,rho,A,I,L,ne,k1,k2)
	#
	# mass and stiffness matrices for a beam using Hermite elements
	# degres of freedom of the element: [y1 rot1 y2 rot2]
	# 
	# E:   elasticity modulus
	# rho: density
	# A:   element section
	# I:   inercia
	# L:   length of the beam
	# ne:  number of elements
	# k1:  added stiffness in vertical direction of first node
	# k2:  added stiffness in vertical direction of last node
	# 
	# javier.cara@upm.es, 
	#

	# local matrices
	# *******************************

	# element length
	le = L/ne

	# matriz de rigidez del elemento
	ke = [12 6*le -12 6*le;
	6*le 4*le^2 -6*le 2*le^2;
	-12 -6*le 12 -6*le;
	6*le 2*le^2 -6*le 4*le^2]
	ke = E*I/(le^3)*ke
		     
	# matriz de masa del elemento (consistente)
	me = [156 22*le 54 -13*le;
	22*le 4*le^2 13*le -3*le^2;
	54 13*le 156 -22*le;
	-13*le -3*le^2 -22*le 4*le^2]
	me = rho*A*le/420*me

	# matrices globales
	# *******************************

	# numero de grados de libertad
	ndof = (ne+1)*2

	K = zeros(ndof,ndof)
	M = zeros(ndof,ndof)
	for e=1:ne
		 dof1 = (e-1)*2 + 1
		 dof2 = (e-1)*2 + 4
		 K[dof1:dof2,dof1:dof2] = K[dof1:dof2,dof1:dof2] + ke
		 M[dof1:dof2,dof1:dof2] = M[dof1:dof2,dof1:dof2] + me
	end

	# muelles en los extremos en direccion y
	K[1,1] = K[1,1] + k1
	K[2*ne+1,2*ne+1] = K[2*ne+1,2*ne+1] + k2

	return K,M

end

# ================================================================
function beam01_test()

	E = 2.1e11
	rho = 7850.
	b = 2.
	h = 0.20
	L = 20.
	ne = 18

	I = 1/12*b*h^3
	A = b*h

	k1 = 100*E*I
	k2 = 100*E*I

	(K,M) = beam01(E,rho,A,I,L,ne,k1,k2)

	# frecuencias naturales del model FEM
	(W,V) = eig(K,M)
	wo = sqrt(W) # (rad/s)

	# frecuencias naturales teoricas
	# wn=Cn*sqrt(E*I/(rho**L^4), donde Cn=(n*pi)^2
	wn = zeros(5)
	for m=1:5
		 wn[m] = (m*pi)^2*sqrt(E*I/(rho*A*L^4))
	end

	return wn,wo[1:5]

end

