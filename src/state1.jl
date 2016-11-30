function state1(m,c,k,f,dt,u0,v0)

	#--------------------------------------------------------------------------
	# Esta funcion calcula la respuesta de un sistema de 1 dgl aplicando el
	# metodo del espacio de los estados
	#
	# Datos de entrada
	# m          masa (kg)
	# k          rigidez (N/m)
	# c          amortiguamiento 
	# f          vector de fuerzas (N)
	# dt         incremento de tiempo
	# u0         posicon inicial
	# v0         velocidad inicial
	#
	# Datos de salida
	# u          desplazamiento de un sistema de 1 gdl sometido a x
	# v          velocidad del sistema
	# a          aceleracion del sistema
	#
	# javier.cara@upm.es, january 2016
	#--------------------------------------------------------------------------

	# numero de instantes de tiempo
	N=length(f)

	# matrices Ac y Bc
	Ac = [0 1;-k/m -c/m]
	Bc = [0;1/m]

	# matrices A y B
	Ad = expm(Ac*dt)
	Bd = (Ad-eye(2))*inv(Ac)*Bc

	# calculo de los estados
	X = zeros(2,N)

	# posicion y velocidad iniciales:
	# en el instante cero tenemos las condiciones iniciales y el primer valor de f
	X[:,1] = [u0;v0]
	for r in 1:N-1
		 X[:,r+1] = Ad*X[:,r] + Bd*f[r]
	end

	u = X[1,:]
	v = X[2,:]
	a = - k/m*u - c/m*v + 1/m*f

	return u,v,a
	
end

# =============================================
function state1_test()

	m = 1.
	w=2*pi # frecuencia propia del sistema
	z=0.02 # tasa de amortiguamiento
	k = m*w^2
	c = 2*z*m*w
	dt = 0.01
	f=randn(1,100)
	u0 = 0.
	v0 = 0.

	u,v,a = state1(m,c,k,f,dt,u0,v0)

	return u,v,a


end
