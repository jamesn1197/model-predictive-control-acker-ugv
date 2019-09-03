#Imports
import numpy as np
import scipy
import scipy.integrate as integrate
from scipy.optimize import minimize
from scipy.linalg import expm
import matplotlib.pyplot as plt

###############################
#...Model Predictive Contol...#
###############################

class MPC:
  def __init__(self,L,h,Np,Q,R):
    self.L=L
    self.h=h
    self.Np=Np
    diagQ=np.kron(np.eye(self.Np),Q)
    diagR=np.kron(np.eye(self.Np),R)
    self.Q=diagQ
    self.R=diagR
  
  def cost(self,du):
    #Something's incompatible here in that the optimizing inputs have to be 1-d and it isn't automatically indexing that way.
    #It keeps randomly not sending the right indexing for the input of the cost function
#     print "Shape of du=",du.shape
#     print "cost du=",du
    self.y=self.rk4(self.x,du.reshape(2,self.Np))[:,1:]
    # print "Shape of y=",self.y.shape
    self.y=self.y.T.reshape(3*self.Np,1)

    #du=du.reshape([self.Np,1])
    #return 0.5*((self.ref-self.W.dot(self.xa)-self.Z.dot(du)).T.dot(self.Q).dot(self.ref-self.W.dot(self.xa)-self.Z.dot(du))+du.T.dot(self.R).dot(du))
    return (self.ref-self.y).T.dot(self.Q).dot(self.ref-self.y)+du.reshape([1,2*self.Np]).dot(self.R).dot(du.reshape([2*self.Np,1]))
  
  def constr0(self,u):
    return u[Np-1]+1
  
  def constr1(self,u):
    return -u[Np-1]+1
  
  def constr2(self,u):
    return u[-1]+0.5
  
  def constr3(self,u):
    return -u[-1]+0.5
  
  def rk4(self,x0,u):
    #4th order Runge-Kutta implementation after collapsing recursive formulation
    
    #x=
    #h=Publishing/Sampling Frequency
    #Np=Number of Predictions
    #u[:,t]=control inputs at the corresponding t
    #h=sample period
    
    u=np.kron(u,np.ones([1,self.Np+1]))
    t_vec=np.linspace(0,self.h*self.Np,self.Np+1)

    #y=np.copy(x)
    x=np.zeros([3,self.Np+1],dtype="float64")
    x[:,0]=x0.reshape(3,)
    #Collapsed RK4 equation
    for t in range(self.Np):
      dx=self.ugvModel(x[:,t],t*self.h,u[:,t])

      #lazy workaround of fencepost error on last iter
      try:
        x[:,t+1]=x[:,t]+(self.h/6*np.array([[2*u[0,t]*np.cos(x[2,t]+self.h*dx[2]/2)*(2+np.cos(self.h*dx[2]/2))],\
                                            [2*u[0,t]*np.sin(x[2,t]+self.h*dx[2]/2)*(2+np.cos(self.h*dx[2]/2))],\
                                            [6*dx[2]]])).reshape(3,)
      except:
        pass
        #print "erk4 fencepost error" if: t==len(tvec)-1

    return x
    
  def optimizer(self,x,R,Q,ref):
    self.x=x.reshape([3,1])
    self.ref=np.kron(np.ones([self.Np,1]),ref)
#     self.xa=xa.reshape([6,1])
#     self.W,self.Z,self.Ca,self.PhiA,self.GammaA=mpc.mpcMats(A,B,C)
#     self.ref=np.kron(np.ones([self.Np,1]),ref)
    bv=[-1,1]
    bgam=[-0.5,0.5]
    bnds=np.kron(np.ones([self.Np,1]),[bv,bgam])
#     print "bnds",bnds
#     print "shape of bounds",bnds.shape
    constr=[{'type':'ineq','fun':self.constr0},\
            {'type':'ineq','fun':self.constr1},\
            {'type':'ineq','fun':self.constr2},\
            {'type':'ineq','fun':self.constr3}]
    init=np.concatenate([.9*np.ones(self.Np,),0.*np.random.rand(self.Np,)])
    print "init=",init
#     init=np.concatenate([0.*np.ones(Np,)])
    dU=minimize(self.cost,init,method='SLSQP',bounds=bnds,)#constraints=constr)
    print(dU.x)
#     print "error",self.ref-self.y
#     print "solved cost:",self.cost(dU.x)
    return dU.x.T
  
  def ugvModel(self,x,t,u):
#     print "x shape",x.shape
#     print "u shape",u.shape
    dx=np.empty(3);
    dx[0]=u[0]*np.cos(x[2]);
    dx[1]=u[0]*np.sin(x[2]);
    dx[2]=u[0]/self.L*np.tan(u[1]);
    
    return dx

########################
#...Operating Script...#
########################
h=0.1
Np=10 #Length of Prediction Horizon
ploss=0. #Probability of lost control signal
closs=0  #Number of consecutive losses (caps at 9)

ref=np.array([-3.,4.,0.]).reshape(3,1)
x0=np.array([3.6,0.,np.pi/2]).T
u0=np.array([1.,0.]).T
tspan=[0,20]
xn=np.zeros([3,int(np.ceil(tspan[1]/h+1))])
xn[0:3,0]=x0
xa=np.zeros([6,int(np.ceil(tspan[1]/h+1))])
u=np.zeros([2,int(np.ceil(tspan[1]/h+1))])#np.kron(np.ones([1,int(np.ceil(tspan[1]/h+1))]),u0.reshape([2,1]))
R=np.diag([0.,0.])
# R=np.diag([0.])
Q=np.diag([4.,4.,0.]) #If I were to impose a cost on theta, it would have to be nonlinear otherwise most initial conditions would throw off the algorithm
xproj=np.zeros([Np,np.ceil(tspan[1]/h+1)],dtype="float64")
yproj=np.zeros_like(xproj)
mpc=MPC(0.2,h,Np,Q,R)

DU=np.zeros([Np,2], dtype='float64')
dU=np.array([0.,0.])

for k in range(1,int(np.ceil(tspan[1]/h+1))):
  print "Loop #:",str(k)
  #Simulate UGVs for iter
  x=integrate.odeint(mpc.ugvModel,xn[:,k-1],np.linspace((k-1)*h,k*h),args=(u[:,k-1],))
  
  #Value at end of time step
  xn[:,k]=x.T[:,-1]  
  
  #C is setup to automatically take the difference of the poses, so r=<0,0,0>
  #r=np.zeros([3,1],dtype="float64")
  event=np.random.rand()
  if event>ploss:
    print "Received"
    DU=mpc.optimizer(xn[:,k],R,Q,ref).reshape(2,Np).T
    print "DU=",DU
    # print "1:",xproj[:,k].shape
    # print "2:",mpc.y[0::3,:].shape
    xproj[:,k-1]=mpc.y[0::3,:].reshape(10,)
    yproj[:,k-1]=mpc.y[1::3,:].reshape(10,)
    # print "prediction=",mpc.y
#     print "solved DU",DU
    closs=0
  else:
    closs+=1
    closs=min(closs,Np-1)
    print "Dropped" 
#   print(DU)
#   dU[0]=np.min([np.max([DU[closs,0],-0.5]),0.5])
#   dU[1]=np.min([np.max([DU[closs,1],-0.1]),0.1])
#   u[0,k]=np.min([np.max([dU[0]+u[0,k-1],-0.5]),1.])
#   u[1,k]=np.min([np.max([dU[1]+u[1,k-1],-0.5]),0.5])
  u[0,k]=DU[closs,0]#+u[1,k-1]
  u[1,k]=DU[closs,1]#+u[1,k-1]
  print "="*30 
  err=xn[0:2,k].reshape([2,1])-ref[0:2,:]
  print err
  print err.T.dot(err)
  if err.T.dot(err)<0.01:
    break
print "Done!",k

#############
#...Plots...#
#############
print(xn.shape)
print(u.shape)

#input plot
plt.plot(u[0,:].T,label="v-mpc")
plt.plot(u[1,:].T,label="gamma-mpc")
plt.xlabel('time')
plt.ylabel('u(t)')
plt.legend()
plt.show()

#xy plot
for k in range(int(np.ceil(tspan[1]/h+1))):
  plt.scatter(xproj[:,k],yproj[:,k],color='red')
  plt.plot(xproj[:,k],yproj[:,k],color='green')
  
plt.plot(xn[0,:],xn[1,:],label="ugv1")
plt.scatter(ref[0],ref[1],color='black')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()

#x plot
# plt.plot(xn[0,:],label="x_ugv1")
# plt.legend()
# plt.show()

#y plot
# plt.plot(xn[1,:],label="y_ugv1")
# plt.legend()
# plt.show()

#theta plot
# plt.plot(xn[2,:],label="th_ugv1")
# plt.legend()
# plt.show()

#error plot
# plt.plot(ref[0]-xn[0,:],label="ex_mpc")
# plt.plot(ref[1]-xn[1,:],label="ey_mpc")
# plt.legend()
# plt.show()

#derivative plot
# plt.plot(xa[6,:],label="dx")
# plt.plot(xa[7,:],label="dy")
# plt.plot(xa[8,:],label="dtheta")
# plt.legend()
# plt.show()

# plt.plot(xn[3,:-1]-xn[3,1:],label="x2mpc")
# plt.plot(xn2[3,:-1]-xn2[3,1:],label="x2pid")
# plt.show()
#print(xn)

