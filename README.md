The purpose of this code is to compute some integral kernels on the sphere S^2 and compare them. They are the following.

* The angular part of the Boltzmann collision kernel corresponding to power-law potentials.
* The kernel of the fractional Laplacian on the sphere.
* The jump kernel of Levy processes subordinate to Browniand motion.

An example of its usage, and the result of some comparison between these kernels can be seen in the Jupyter notebooks that are provided. In [this notebook](https://github.com/luissilvestre/collisionkernel/blob/main/fractional_laplacian_test.ipynb) we compare the collision kernel with the spherical fractional Laplacian. In [this notebook](https://github.com/luissilvestre/collisionkernel/blob/main/subordinate_test.ipynb) we compare the collision kernel with the kernel of a subordinate process.

The code is written using the Julia programming languaga. If you have no idea how to make this work, you may still run the code with this link [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/luissilvestre/collisionkernel/main?labpath=subordinate_experiments.ipynb). Beware that the computations on the cloud server are very slow.
