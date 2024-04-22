import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scp
def boolean_gate(weights, r):
    num_weights = len(weights)
    min_weight = min(weights)
    max_weight = max(weights)
    boolean_value = r * min_weight + (1 - r) * max_weight
    return boolean_value
simulation_time = 1500  # Define the simulation time
dt = 0.1
Proteins = ["cGMP", "PKG", "CREB", "BAD", "VASP"]
Number_of_Protein = len(Proteins)  # Define the number of proteins
activities = np.zeros((simulation_time, Number_of_Protein))
activities[0, :] = [1, 0.1, 0.1, 0.1, 0.1]
protein = np.zeros(Number_of_Protein)
cGMPweights = [1] #3rd one for PKG
 #3, 4, 5 are for CREB, BAD and VASP respectively.
sigma = np.full(5, 0)

downstream = np.zeros(Number_of_Protein)
upstream = np.zeros(Number_of_Protein)
new_protein = np.zeros(Number_of_Protein)
r = 0.75

def diff_Inhi_Proli(PKGweights):
    for k in range(simulation_time-1):
        for i in range(Number_of_Protein):

            if activities[k, i]<=0:
                    continue
            else:
                if i==0:

                    wiener = np.random.normal(0, dt)

                    activities[k+1, i] = activities[k, i]-(1*boolean_gate(cGMPweights, 0.75)*activities[k, i])*dt+sigma[i]*activities[k, i]*wiener
                    if activities[k+1, i]<=0: activities[k+1, i] = 0

                elif i==1:

                    wiener = np.random.normal(0, dt)

                    activities[k+1, i]=activities[k, i]+(1*boolean_gate([cGMPweights[0]], 0.75)*activities[k, 0]-1*boolean_gate(PKGweights, 0.75)*activities[k, i])*dt+sigma[i]*activities[k, i]*wiener
                    if activities[k+1, i]<=0: activities[k+1, i] = 0
                else:

                    #print(activities[k, i])
                    wiener = np.random.normal(0, dt) 
                    activities[k+1, i]=activities[k, i]+(1*boolean_gate([PKGweights[i-2]], 0.75)*activities[k, 1])*dt+sigma[i]*activities[k, i]*wiener
                    if activities[k+1, i]<=0: activities[k+1, i] = 0
    Inhibition = activities[:, 3]*activities[:, 2]
    return(abs(Inhibition[simulation_time-1] - activities[simulation_time-1, 4]))

    
bounds = [(0, 1),(0, 1),(0, 1)]
pkg_initial = [0.5, 0.5, 0.5]
# Constraint function
#def constraint(p):
 #   return p[0]+p[1]+p[2]+p[3]-1 

#cons = ({'type': 'eq', 'fun': constraint})

results = scp.minimize(diff_Inhi_Proli, pkg_initial, bounds=bounds, method='SLSQP')
print(results.x, results.fun, results.x[0]*results.x[1]/results.x[2])
              
Inhibition = activities[:, 3]*activities[:, 2]
plt.xlabel("Generations")
plt.ylabel("Activity")
plt.plot(activities[:, 0], "b", label = "cGMP")
plt.plot(activities[:, 1], "c", label = "PKG")
plt.plot(activities[:, 2], "r", label = "CREB")
plt.plot(activities[:, 3], "y", label = "BAD")
plt.plot(activities[:, 4], "g", label = "VASP")
plt.title("Optimised Evolution")
plt.legend()
plt.savefig("Optimised Evolution.png", dpi=300)
plt.close()

plt.xlabel("Generations")
plt.ylabel("Activity")
plt.plot(Inhibition, "black", label = "Inhibition of Apoptosis")
plt.plot(activities[:, 4], "g", label = "Cell Proliferation")
plt.title("Optimised Inhibition and Cell Proliferation")
plt.legend()
plt.savefig("Optimised Inhibition and Cell Proliferation.png", dpi=300)
plt.close()