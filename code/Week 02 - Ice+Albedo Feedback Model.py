import numpy as np
import matplotlib.pyplot as plt

# ask user for inputs
L, albedo, n_iterations = input("").split()
L, albedo, n_iterations = [ float(L), float(albedo), int(n_iterations) ]

L_range = (1200, 1600)
#n_iterations = 100
# albedo = 0.65 
# slope (m) and intersection (b) of lines for albedo/ice latitude vs temperature
m_albedo = -1e-2
b_albedo = 2.8
m_ice_lat = 1.5
b_ice_lat = -322.5
max_albedo = 0.65
min_albedo = 0.15

sigma = 5.670374419e-8 # stefan-boltzmann constant

x = []
y = []

plot_type = 'iter_down' # 'L', 'iter_up', 'iter_down'


#'iter up'
#L = L_range[0] #[L_range[1]]
while L < L_range[1]+1:
    for i in range(n_iterations):
        T = (L*(1-albedo)/(4*sigma))**(1/4)
        albedo = T*m_albedo + b_albedo
        albedo = min(albedo, max_albedo)
        albedo = max(albedo, min_albedo)
        #ice_lat = T*m_ice_lat + b_ice_lat
        #ice_lat = min(ice_lat, 90)
        #ice_lat = max(ice_lat, 0)
        
        if plot_type  == 'iter_up': #is 'iter_up':
            x.append(i)
            y.append(T)
    if plot_type == 'iter_up': #is 'iter_up':
        x.append(np.nan)
        y.append(np.nan)
    if plot_type  == 'L': #is 'L':
        x.append(L)
        y.append(T)
    L += 10

# iter_down
#L = L_range[1]
while L > L_range[0]-1:
    for i in range(n_iterations):
        T = (L*(1-albedo)/(4*sigma))**(1/4)
        albedo = T*m_albedo + b_albedo
        albedo = min(albedo, max_albedo)
        albedo = max(albedo, min_albedo)
        #ice_lat = T*m_ice_lat + b_ice_lat
        #ice_lat = min(ice_lat, 90)
        #ice_lat = max(ice_lat, 0)
        
        if plot_type == 'iter_down': #is 'iter_down':
            x.append(i)
            y.append(T)
    if plot_type == 'iter_down': #is 'iter_down':
        x.append(np.nan)
        y.append(np.nan)
    if plot_type == 'L': #is 'L':
        x.append(L)
        y.append(T)
    L -= 10
    
plt.plot(x,y)
plt.xlabel(plot_type)
plt.ylabel('T')
plt.title('T as a function of L (and the albedo)')
plt.show()

print(T, albedo)