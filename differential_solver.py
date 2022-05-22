import numpy as np


#Eulers - explicit method
def euler(f,x0,t,h):
    x = np.zeros(len(t))
    x[0] = x0

    for i in range(0,len(t)-1):
        x[i+1] = x[i] +h*f(t[i],x[i])
    
    return x

def implicit_meth(f,x0,t,h):
    x = np.zeros(len(t))
    x[0] = x0

    for i in range(0,len(t)-1):
        x[i+1] = x[i] +h*f(t[i+1],x[i+1])
    
    return x


def predictor(f,y_0,x_0,n,h):
	y=[]
	x=[]
	yc=[]
	for j in range(0,n+1):
		y.append(0)
		x.append(0)
		yc.append(0)
	y[0]=y_0
	x[0]=x_0
	

	for i in range(0,n):
		x[i]=x[0]+i*h
		k1 = h*f(y[i],x[i])
		k2 = h*f(y[i+1],x[i]+h)
		yc[i+h] = y[i]+(k1+k2)/2

	return yc	

#RK$ 
def RungKutta4(f,u0,xMin,xMax,n):
    #where f is dy/dx function
    #u0 is array of initial conditions based on copupled equations in f
    #t is range of arrays with step size of dt 

    x = np.linspace(xMin,xMax,n)
    nx = len(x)

    nu = len(u0)
    u = np.zeros((nx,nu))
    u[0] = u0

    for k in range(nx-1):
        dx = (xMax-xMin)/n
        k1 = dx*f(u[k], x[k])
        k2 = dx*f(u[k]+(k1/2), x[k]+(dx/2))
        k3 = dx*f(u[k]+(k2/2),x[k]+(dx/2))
        k4 = dx*f( u[k]+k3,x[k]+dx)

        du = (k1 + 2*k2 + 2*k3 +k4)/6

        u[k+1] =u[k]+du 

    #return the increment of value with each step as an array     
    return u


#RK4 for shooting method 

def  RK4(f,g,x0,y0,z0,xn,h):
	#creating arrays to store large date for plot
	x=[]
	y=[]
	xl = -xn
	#checking with the boundary condition adn thus applying condition
	while x0>xl:
		k1 = -h*f(x0,y0,z0)
		#if another equation 
		l1 = -h*g(x0,y0,z0)

		k2 = -h*f(x0-h/2,y0+k1/2,z0+l1/2)
		l2 = -h*g(x0-h/2,y0+k1/2,z0+l1/2)

		k3 = -h*f(x0-h/2,y0+k2/2,z0+l2/2)
		l3 = -h*g(x0-h/2,y0+k2/2,z0+l2/2)

		k4 = -h*f(x0-h/2,y0+k3/2,z0+l3/2)
		l4 = -h*g(x0-h/2,y0+k3/2,z0+l3/2)

		#x0,y0,z0 are boundary conditions

		z0 = z0+(l1 + 2*(l2+l3)+l4)/6
		y0 = y0+(k1 + 2*(k2+k3)+k4)/6
		x0 = x0 -h
		#append it to the array
		y.append(y0)
		x.append(x0)

	while x0<xn:
		k1 = h*f(x0,y0,z0)
		#if another equation 
		l1 = h*g(x0,y0,z0)

		k2 = h*f(x0+h/2,y0+k1/2,z0+l1/2)
		l2 = h*g(x0+h/2,y0+k1/2,z0+l1/2)

		k3 = h*f(x0+h/2,y0+k2/2,z0+l2/2)
		l3 = h*g(x0+h/2,y0+k2/2,z0+l2/2)

		k4 = h*f(x0+h/2,y0+k3/2,z0+l3/2)
		l4 = h*g(x0+h/2,y0+k3/2,z0+l3/2)

		z0 = z0+(l1 + 2*(l2+l3)+l4)/6
		y0 = y0+(k1 + 2*(k2+k3)+k4)/6
		x0 = x0 +h
		y.append(y0)
		x.append(x0)
	N = len(x)	
	for i in range(0,N):	
		with open('data.txt','a') as l:
			print(x[i],',',y[i],file=l)		

	return x,y,y0
	
def shooting_method(f,g,x0,y0,z0,z1,xn,yn,h):
	#solving differential eq using RK4 method for mentioned boundary conditions
	a,b,c = RK4(f,g,x0,y0,z0,xn,h)
	p,q,r = RK4(f,g,x0,y0,z1,xn,h)
	
	if (c-yn)>0.001 and (yn-r)>0.001:
		xl=(z0-z1)*(yn-r)
		z= z1 + xl/(c-r)
		s,t,u = RK4(f,g,x0,y0,z,xn,h)

		if (u-yn)>0.001:
			#using Lagrange's interpolation formula
			z2=z1+(z-z1)*(yn-r)/(u-r)
			shooting_method(f,g,x0,y0,z2,z1,xn,yn,h)
		if (yn-u)>0.001:
			z2=z+(z0-z)*(yn-u)/(c-u)
			shooting_method(f,g,x0,y0,z2,z1,xn,yn,h)
		if abs(c-yn)<0.001:
			return s,t,u,z
	if abs(c-yn)<0.001 :
		return a,b,c,z0
	if abs(yn-r)<0.001 :
		return p,q,r,z1
	return z	


