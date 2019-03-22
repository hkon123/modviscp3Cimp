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
        self.imei, = plt.plot(self.t,self.en)
        plt.ylim(-600,0)
        plt.xlim(0,100000)
        anim = FuncAnimation(fig, self.update, frames = 100, repeat = False, interval = 1, blit = True)
        plt.show()

    def update(self, i):
        self.count+=self.sweeps
        if i ==0:
            self.phi, ensum = example_test.f(np.random.uniform(-0.01,0.01,size = [100,100]),self.sweeps)
            self.en = np.append(self.en,ensum)
        else :
            self.phi, ensum = example_test.f(self.phi,self.sweeps)
            self.en = np.append(self.en,ensum)
        self.t = np.append(self.t,self.count)
        self.imei.set_data(self.t,self.en)
        return [self.imei]


A=Anim(100)
