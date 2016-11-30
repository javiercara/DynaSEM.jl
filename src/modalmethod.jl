function modalmethod(M,K,zm,F,nm,dt)

	# eigenvalues and eigenectors
	D,V = eig(K,M)

	# modal matrices
	Mm = V'*M*V
	Km = V'*K*V
	Fm = V'*F

	ndof,nt = size(F)
	
	qm=zeros(ndof,nt)
	q1m=zeros(ndof,nt)
	q2m=zeros(ndof,nt)
	for modo in 1:nm
		  m = Mm[modo,modo]
		  k = Km[modo,modo]
		  c = 2*sqrt(m*k)*zm[modo]
		  f = Fm[modo,:]
		  
		  ym,vm,am = state1(m,c,k,f,dt,0.,0.)
		  qm[modo,:] = ym
		  q1m[modo,:] = vm
		  q2m[modo,:] = am
	end

	# deshacemos el cambio
	q = V*qm
	q1 = V*q1m
	q2 = V*q2m
	
	return q,q1,q2
	
end

# ========================================================
function modalmethod_test()
	include("FEMmodel.jl")
	E = 2.1e11
	rho = 7850.
	b = 2.
	h = 0.20
	L = 20.
	ne = 18

	I=1/12*b*h^3
	A=b*h

	k1=100*E*I
	k2=100*E*I

	(K,M) = FEMmodel(E,rho,A,I,L,ne,k1,k2)

	ndof = size(K,1)
	zm = 0.02*ones(ndof) # modal dampings
	F = randn(ndof,100)
	
	nm=3
	
	dt = 0.1
	
	q,q1,q2 = modalmethod(M,K,zm,F,nm,dt)
	
	return q,q1,q2
	
end
