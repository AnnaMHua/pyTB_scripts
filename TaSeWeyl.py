from pythtb import * # import TB model class
import matplotlib.pyplot as plt

# 3D
dimr = 3
dimk = 3

lat=[[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]
orb=[[0, 0, 0], [0, 0, 0.5], [0.5, 0.5, 0], [0.5, 0.5, 0.5]]

my_model=tb_model(dim_k=dimk,dim_r=dimr,lat=lat,orb=orb)

# set model parameters
v=1
u_1 = 0.1
u_1t = 0.2
u_2 = 0.1

# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
# H1 term: intra-chain hopping
#   hopping along z-axis within single chain
#       I-chain intra chain
my_model.set_hop(-v, 0, 1, [ 0, 0, 0])
my_model.set_hop(-v, 0, 1, [ 0, 0, -1])
#       II-chain intra chain
my_model.set_hop(-v, 2, 3, [ 0, 0, 0])
my_model.set_hop(-v, 2, 3, [ 0, 0, -1])

# H2 term in paper: inter-chain hopping
# hopping within the sublattice
#   I-chain sublattice hopping
#   x-direction
my_model.set_hop(-u_1, 0, 1, [ 1, 0, 0])
my_model.set_hop(-u_1, 0, 1, [ 1, 0, -1])
my_model.set_hop(-u_1, 0, 1, [ -1, 0, 0])
my_model.set_hop(-u_1, 0, 1, [ -1, 0, -1])

#   y-direction
my_model.set_hop(-u_1t, 0, 1, [0, 1, 0])
my_model.set_hop(-u_1t, 0, 1, [0, 1, -1])
my_model.set_hop(-u_1t, 0, 1, [0, -1,  0])
my_model.set_hop(-u_1t, 0, 1, [0, -1,  -1])

# II-chain sublattice hopping
#   x-direction
my_model.set_hop(-u_1t, 2, 3, [ 1, 0, 0])
my_model.set_hop(-u_1t, 2, 3, [ 1, 0, -1])
my_model.set_hop(-u_1t, 2, 3, [ -1, 0, 0])
my_model.set_hop(-u_1t, 2, 3, [ -1, 0, -1])

#   y-direction
my_model.set_hop(-u_1, 2, 3, [0, 1, 0])
my_model.set_hop(-u_1, 2, 3, [0, 1, -1])
my_model.set_hop(-u_1, 2, 3, [0, -1,  0])
my_model.set_hop(-u_1, 2, 3, [0, -1,  -1])

# H3 term in paper
my_model.set_hop(-u_2, 0, 0, [1, 1, 0])
my_model.set_hop(u_2, 0, 0, [1, -1, 0])
my_model.set_hop(u_2, 1, 1, [1, 1, 0])
my_model.set_hop(-u_2, 1, 1, [1, -1, 0])

my_model.set_hop(u_2, 2, 2, [1, 1, 0])
my_model.set_hop(-u_2, 2, 2, [1, -1, 0])
my_model.set_hop(-u_2, 3, 3, [1, 1, 0])
my_model.set_hop(u_2, 3, 3, [1, -1, 0])

# define a path in k-space

# Gamma, Z, A, M, Gamma, Z, R, A, Z
path=[[0., 0., 0.], [0., 0., 0.5], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0, 0, 0],
       [0, 0, 0.5], [0, 0.5, 0.5], [0.5,0.5,0.5], [0, 0, 0.5],]
k_label=[r"$\Gamma$",r"$Z$", r"$A$", r"$M$", r"$\Gamma$",r"$Z$", r"R", r"A", r"Z"]

# # Z, R, A
# path=[[0, 0, 0.5], [0, 0.5, 0.5], [0.5,0.5,0.5]]
# k_label=[r"$Z$", r"R", r"A"]

(k_vec,k_dist,k_node)=my_model.k_path(path,5000)

# solve model
evals=my_model.solve_all(k_vec)

# plot band structure
fig, ax = plt.subplots()
ax.plot(k_dist,evals[0])
ax.plot(k_dist,evals[1])
ax.plot(k_dist,evals[2])
ax.plot(k_dist,evals[3])
ax.set_title("TaSe Weyl model band structure")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")
ax.set_xticks(k_node)
ax.set_xticklabels(k_label)
ax.set_xlim(k_node[0],k_node[-1])
for n in range(len(k_node)):
  ax.axvline(x=k_node[n], linewidth=0.5, color='k')

ax.text(0.7, 1.5, f"v={v}\nu_1={u_1}\nu_1t={u_1t}\nu_2={u_2}")
fig.tight_layout()
fig.savefig("TaSeWeyl.pdf")