function chirp_lin(f0,t1,f1,fs,phase)
	#
	# f0: frequency at t=0 (Hz)
	# f1: frequency at t=t1 (seconds)
	# fs: sampling frequency (Hz)
	# phase: phase at t=0 (degrees)
	#
	# variable frequency - linear
	# f(t) = f0 + b*t, where b = (f1-f0)/t1
	# 
	# javier.cara@upm.es, 2016 - january

	nt = Int(t1*fs)
	t = (0:nt-1)*(1/fs)

	b = (f1-f0)/t1
	x = Array(Float64,1,nt)
	for k in 1:nt
		  x[k] = sin(2*pi*(f0*t[k] + b/2*t[k]^2) + phase/360*2*pi)
	end

	return x

end

# ***************************************************************
function chirp_lin_test()
	
	f0=0.
	f1=20.
	t1=50.
	fs=100.
	phase=90.

	nt = Int(t1*fs)
	t = (0:nt-1)*(1/fs)
		
	x = chirp_lin(f0,t1,f1,fs,phase)
		
	return x,t
end

# **************************************************************
function chirp_log(f0,t1,f1,fs,phase)
	#
	# f0: frequency at t=0 (Hz)
	# f1: frequency at t=t1 (seconds)
	# fs: sampling frequency (Hz)
	# phase: phase at t=0 (degrees)
	#
	# variable frequency - linear
	# f(t) = f0*b^t, where b = (f1/f0)(1/t1)
	# 
	# javier.cara@upm.es, 2016 - january
	
	if f0 == 0
	# f0 can't be zero
		f0 = 1.0e-6
	end
	
	nt = Int(t1*fs)
	t = (0:nt-1)*(1/fs)

	b = (f1/f0)^(1/t1)
	x = Array(Float64,1,nt)
	for k in 1:nt
		  x[k] = sin(2*pi*f0*(b^t[k]-1)/log(b) + phase/360*2*pi)
	end

	return x

end

# ***************************************************************
function chirp_log_test()
	
	f0=0.
	f1=20.
	t1=50.
	fs=100.
	phase=90.

	nt = Int(t1*fs)
	t = (0:nt-1)*(1/fs)
		
	x = chirp_log(f0,t1,f1,fs,phase)
		
	return x,t
end

