import numpy as npclass World:        def __init__(self, change=None, name=''):        self.name = name        self.change = change        self.time = []        self.dt = None        self.pCO2 = []        self.RF = []        self.RF_CO2 = []        self.RF_masking = []        self.Teq = []        self.T = []    def initialize_time(self, timespan=[1900, 2100], nsteps=200):        self.dt = (timespan[1] - timespan[0])/nsteps        self.time = list(np.arange(timespan[0], timespan[1]+self.dt, self.dt))            def compute_CO2(self, pCO2_ini=290, A=0.0225):        # initialize        self.pCO2.append(pCO2_ini)        # integrate        for t in self.time[1:]:            if self.change is None or t < self.change:                self.pCO2.append(                    280 + (1+A*self.dt)*(self.pCO2[-1]-280))            else:                self.pCO2.append(                    self.pCO2[-1] + (0.01*self.dt)*(340-self.pCO2[-1]))                    def compute_RF(self, aerosol_Wm2_now=-0.75, C=4):        # calculate B        B = aerosol_Wm2_now*self.dt/(self.pCO2[self.time.index(2014)]                                     - self.pCO2[self.time.index(2013)])        masking2015 = -float('inf')        for i, t in enumerate(self.time[1:]):            self.RF_CO2.append(C*np.log(self.pCO2[i+1]/280)/np.log(2))            if self.change is None or t < self.change:                self.RF_masking.append(                    max(B*(self.pCO2[i+1]-self.pCO2[i])/self.dt, masking2015))                if t == 2015:                    masking2015 = self.RF_masking[-1]            else:                self.RF_masking.append(0)            self.RF.append(self.RF_CO2[-1]+self.RF_masking[-1])        # add 0th values to get same length as time        self.RF_CO2 = [self.RF_CO2[0]] + self.RF_CO2        self.RF_masking = [self.RF_masking[0]] + self.RF_masking        self.RF = [self.RF[0]] + self.RF     def compute_T(self, C=4, climate_sensitivity_2x=3):        self.Teq.append(self.RF[0]/C*climate_sensitivity_2x)        self.T.append(self.RF[0]/C*climate_sensitivity_2x)        for i, t in enumerate(self.time[1:]):            self.T.append(self.T[-1]+(self.Teq[-1]-self.T[-1])/(20/self.dt))            self.Teq.append(self.RF[i+1]/C*climate_sensitivity_2x)                def run(self):        self.initialize_time()        self.compute_CO2()        self.compute_RF()        self.compute_T()            def axes(self, ax):        ax[2].set_xlabel('years')        ax[0].set_ylabel('pCO2 [ppm]')        ax[0].grid()        ax[1].set_ylabel('RF [W/m2]')        ax[1].grid()        ax[2].set_ylabel('T [C]')        ax[2].grid()        ax[0].set_xlim([min(self.time), max(self.time)])                    def plot(self, ax, color='r'):        ax[0].plot(self.time, self.pCO2, color=color,                   label=self.name + ': pCO2')        ax[1].plot(self.time, self.RF_CO2, color=color, ls='--',                   label=self.name + ': RF_CO2')        ax[1].plot(self.time, self.RF_masking, ls=':', color=color,                    label=self.name + ': RF_masking')        ax[1].plot(self.time, self.RF, ls='-', color=color,                    label=self.name + ': RF_total')         ax[2].plot(self.time, self.T, ls='-', color=color,                    label=self.name + ': T')        ax[2].plot(self.time, self.Teq, ls='--', color=color,                    label=self.name + ': Teq')            def legends(self, ax):        ax[0].legend(prop={'size': 6})        ax[1].legend(prop={'size': 6})        ax[2].legend(prop={'size': 6})         if __name__ == '__main__':    plot = True    # scenarios    world_bau = World(change=None, name='business as usual')    world_bau.run()    world_wou = World(change=2016, name='without us after 2016')    world_wou.run()    #plots    if plot is True:        import matplotlib.pyplot as plt        fig, ax = plt.subplots(3,1,sharex=True, figsize=(8,8))        world_bau.axes(ax)        world_bau.plot(ax, color='b')        world_wou.plot(ax, color='r')            world_bau.legends(ax)        plt.show()