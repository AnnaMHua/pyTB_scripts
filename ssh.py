from pythtb import * # import TB model class
import matplotlib.pyplot as plt

lat=[[1.0]]
orb=[[0.0],[0.5]]

my_model=tb_model(1,1,lat,orb)

t = 1.0
delta = 0.1

my_model.set_hop(-(t+delta), 0, 1, [1.0])
my_model.set_hop(-(t-delta), 0, 1, [0])

# define a path in k-space
(k_vec,k_dist,k_node)=my_model.k_path('fullc',1000)
k_label=[r"$0$",r"$\pi$", r"$2\pi$"]

# solve model
evals=my_model.solve_all(k_vec)

# plot band structure
fig, ax = plt.subplots()
ax.plot(k_dist,evals[0])
ax.plot(k_dist,evals[1])
ax.set_title("SSH model band structure")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")
ax.set_xticks(k_node)
ax.set_xticklabels(k_label)
ax.set_xlim(k_node[0],k_node[-1])
for n in range(len(k_node)):
  ax.axvline(x=k_node[n], linewidth=0.5, color='k')
fig.tight_layout()
fig.savefig("SSH_band.pdf")