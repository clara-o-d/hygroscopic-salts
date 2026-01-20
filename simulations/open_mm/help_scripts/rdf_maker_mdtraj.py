import matplotlib.pyplot as plt
import mdtraj
import numpy as np

traj = mdtraj.load("trajectory.dcd", top="start_Mg_Cl_fixed.pdb")

# pairs = np.array([[1, 2]])  # atom indices

pairs = traj.topology.select_pairs("name O", "name O")


r, gr = mdtraj.compute_rdf(
    traj,
    pairs=pairs,
    r_range=(0.1, 1.5),
    bin_width=0.002,
)

plt.close("rdf")
fig, ax = plt.subplots(num="rdf")
ax.plot(r, gr)
ax.set_xlabel("Radius [nm]")
ax.set_ylabel("O-O Radial distribution function g(r)")
fig.savefig("rdfOs.png", dpi=300, bbox_inches="tight")
