import example_test
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation




class Anim(object):




    def chaan_hillard(self,sweeps):
        self.sweeps = sweeps
        self.count = 0
        fig, ax = plt.subplots()
        self.en = np.array((0))
        self.t = np.array((0))
        self.energy = np.zeros((100,100))
        self.phi = np.random.uniform(-0.01,0.01,size = [100,100])
        self.initializePhi()
        plt.subplot(1,3,1)
        self.im=plt.imshow(self.phi, cmap='plasma',vmin = -1, vmax = 1)
        plt.colorbar()
        plt.subplot(1,3,2)
        self.energyim=plt.imshow(self.energy, cmap='plasma')
        plt.colorbar()
        plt.subplot(1,3,3)
        self.imei, = plt.plot(self.t,self.en)
        plt.ylim(-600,0)
        plt.xlim(0,100000)
        anim = FuncAnimation(fig, self.update, frames = 1000, repeat = False, interval = 1, blit = True)
        plt.show()

    def update(self, i):
        self.count+=self.sweeps
        self.phi,self.energy, ensum = example_test.f(self.phi,self.sweeps)
        self.en = np.append(self.en,ensum)
        self.t = np.append(self.t,self.count)
        self.imei.set_data(self.t,self.en)
        self.im.set_array(self.phi)
        self.energyim.set_array(self.energy)
        return [self.im,self.energyim,self.imei]


    def initializePhi(self):
        for i in range(100):
            for j in range(100):
                if np.random.uniform() >=0.5:
                    self.phi[i,j] = 0.5#+np.random.uniform(-0.01,0.01)
                else:
                    self.phi[i,j] = -0.5#+np.random.uniform(-0.01,0.01)

    def jacobi(self, treshold, rho):
        phi,conv,Exj,Eyj,value = example_test.j(rho, treshold)
        print(value)
        y = np.zeros(value)
        for i in range(value):
            y[i] = conv[i]
        cut = np.zeros((50,50))
        xcutj = np.zeros((50,50))
        ycutj = np.zeros((50,50))
        for i in range(50):
            for j in range(50):
                    cut[i,j] =  phi[i,j,25]
                    xcutj[i,j] =  Exj[i,j,1]
                    ycutj[i,j] =  Eyj[i,j,1]
        plt.quiver(xcutj,ycutj)#, scale_units='xy')
        plt.show()
        im = plt.imshow(cut, cmap = 'plasma')
        plt.colorbar()
        plt.show()
        plt.plot(np.arange(0,value,1),y)
        plt.show()

    def gauss_seidel(self,treshold,rho):
        phig,convg,Exg,Eyg,valueg = example_test.g(rho,treshold)
        print(valueg)
        yg = np.zeros(valueg)
        for i in range(valueg):
            yg[i] = convg[i]
        cutg = np.zeros((50,50))
        xcutg = np.zeros((50,50))
        ycutg = np.zeros((50,50))
        for i in range(50):
            for j in range(50):
                    cutg[i,j] =  phig[i,j,25]
                    xcutg[i,j] =  Exg[i,j,1]
                    ycutg[i,j] =  Eyg[i,j,1]
        plt.quiver(xcutg,ycutg)#, scale_units='xy')
        plt.show()
        im = plt.imshow(cutg, cmap = 'plasma')
        plt.colorbar()
        plt.show()
        plt.plot(np.arange(0,valueg,1),yg)
        plt.show()


rho = np.zeros((50,50,50))
rho[25,25,25] = 1
A=Anim().jacobi(0.001, rho)
