import example_test
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation




class Anim(object):

    def __init__(self, sweeps):
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


A=Anim(100)
