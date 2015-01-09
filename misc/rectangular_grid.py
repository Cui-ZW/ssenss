import h5py
import numpy as np

def new_parameters():
	param_type = np.dtype([("n", "i"), ("nel", "i")])
	return np.zeros((1,), dtype = param_type)

def new_elements(n, nel):
	element_type = np.dtype([("Neighbours", "i", (6,)),
			("Neighbour rotation", "i", (6,)),
			("Corner size", "i", (8,)),
			("Corner neighbours", "i", (8,)),
			("Edge size", "i", (12,)),
			("Edge neighbours", "i", (12,)),
			("X", "double", (n**3,)),
			("Y", "double", (n**3,)),
			("Z", "double", (n**3,))])
	return np.zeros((nel,), dtype = element_type)

def write_dataset(fp, dat, name):
	fp.create_dataset(name, data = dat)

def write_parameters(fp, n, nel):
	dat = new_parameters()
	dat[0]["n"] = n
	dat[0]["nel"] = nel
	write_dataset(fp, dat, "Parameters")

def write_elements(fp, n, nel, grid_func, ind_order, xmin, xmax, bc):
	dat = new_elements(n, nel[0] * nel[1] * nel[2])
	set_element_values(dat, nel, grid_func, ind_order, xmin, xmax)
	set_boundary_conditions(dat, nel, ind_order, bc)
	write_dataset(fp, dat, "Element array")

def set_element_values(dat, nel, grid_func, ind_order, xmin, xmax):
	for k in range(nel[2]):
		for j in range(nel[1]):
			for i in range(nel[0]):
				ind = np.array([i, j, k])
				el = elnum(ind[ind_order], nel[ind_order])
				(x, y, z) = gridpoints(ind, nel, xmin, xmax,
						grid_func)
				dat[el]["X"][:] = x
				dat[el]["Y"][:] = y
				dat[el]["Z"][:] = z
				dat[el]["Neighbours"][:] = neighbours(ind, nel,
						ind_order)
				dat[el]["Neighbour rotation"][:] = np.zeros(6)
				dat[el]["Edge neighbours"][:] = edge_neighbours(
						ind, nel, ind_order)
				dat[el]["Edge size"][:] = np.ones(12)
				dat[el]["Corner neighbours"][:] = corner_neighbours(
						ind, nel, ind_order)
				dat[el]["Corner size"][:] = np.ones(8)

def set_boundary_conditions(dat, nel, ind_order, bc):
	for s in range(6):
		if bc[s] < 0:
			k = s // 2
			l = s % 2
			for j in range(np.roll(nel, -k)[2]):
				for i in range(np.roll(nel, -k)[1]):
					ind = np.roll(np.array([l *
							(nel[k] - 1), i, j]), k)
					el = elnum(ind[ind_order],
							nel[ind_order])
					set_boundary_type(dat, el, s, bc[s])

def set_boundary_type(dat, el, surf, bc):
	edge_ind = [[0, 2, 4, 6], [1, 3, 5, 7], [0, 1, 8, 10], [2, 3, 9, 11],
			[4, 5, 8, 9], [6, 7, 10, 11]]
	corner_ind = [[0, 2, 4, 6], [1, 3, 5, 7], [0, 1, 4, 5], [2, 3, 6, 7],
			[0, 1, 2, 3], [4, 5, 6, 7]]
	dat[el]["Neighbours"][surf] = bc
	dat[el]["Edge neighbours"][edge_ind[surf]] = bc
	dat[el]["Edge size"][edge_ind[surf]] = 0
	dat[el]["Corner neighbours"][corner_ind[surf]] = bc
	dat[el]["Corner size"][corner_ind[surf]] = 0

def elnum(i, n):
	return i[0] + i[1] * n[0] + i[2] * n[0] * n[1]

def perind(i, n):
	return i - (i // n) * n

def gridpoints(i, n, xmin, xmax, grid_func):
	(xa, xb) = grid_func[0](i[0], n[0], xmin[0], xmax[0])
	(ya, yb) = grid_func[1](i[1], n[1], xmin[1], xmax[1])
	(za, zb) = grid_func[2](i[2], n[2], xmin[2], xmax[2])
	xgrid = xa + (xb - xa) * np.array([0, 1, 0, 1, 0, 1, 0, 1])
	ygrid = ya + (yb - ya) * np.array([0, 0, 1, 1, 0, 0, 1, 1])
	zgrid = za + (zb - za) * np.array([0, 0, 0, 0, 1, 1, 1, 1])
	return (xgrid, ygrid, zgrid)

def equidistant_grid(i, n, xmin, xmax):
	return (xmin + i * (xmax - xmin) / n,
			xmin + (i + 1) * (xmax - xmin) / n)

def chebyshev_node(i, n, xmin, xmax):
	tmp = np.cos((2.0 * (i + 1) - 1.0) * np.pi / (2.0 * n))
	return 0.5 * (xmin + xmax) - 0.5 * (xmax - xmin) * tmp

def chebyshev_grid(i, n, xmin, xmax):
	cn = [chebyshev_node(j, n, xmin, xmax) for j in [i - 1, i, i + 1]]
	a = 0.5 * (cn[0] + cn[1])
	b = 0.5 * (cn[1] + cn[2])
	if i == 0:
		a = xmin
	elif i == (n - 1):
		b = xmax
	return (a, b)

def neighbours(i, n, ind_order):
	neigh = np.zeros(6)
	for j in range(6):
		di = np.roll(np.array([2 * (j % 2) - 1, 0, 0]), j // 2)
		neigh[j] = elnum(perind(i + di,	n)[ind_order], n)
	return neigh

def edge_neighbours(i, n, ind_order):
	neigh = np.zeros(12)
	for j in range(12):
		di1 = np.roll(np.array([2 * (j % 2) - 1, 0, 0]), j // 8)
		di2 = np.roll(np.array([2 * ((j // 2) % 2) - 1, 0, 0]),
				2 - (11 - j) // 8)
		neigh[j] = elnum(perind(i + di1 + di2, n)[ind_order], n)
	return neigh

def corner_neighbours(i, n, ind_order):
	neigh = np.zeros(8)
	for j in range(8):
		di = np.array([2 * (j % 2) - 1, 2 * ((j // 2) % 2) - 1,
				2 * ((j // 4) % 2) - 1])
		neigh[j] = elnum(perind(i + di, n)[ind_order], n)
	return neigh

def read_user_input():
	print("Rectangular grid generator for SSENSS")
	print("")
	print("Please specify grid parameters")

	xmin = float(input("Minimum x value, xmin = "))
	xmax = float(input("Maximum x value, xmax = "))
	ymin = float(input("Minimum y value, ymin = "))
	ymax = float(input("Maximum y value, ymax = "))
	zmin = float(input("Minimum z value, zmin = "))
	zmax = float(input("Maximum z value, zmax = "))
	nelx = int(input("Number of elements in x direction, nelx = "))
	nely = int(input("Number of elements in y direction, nely = "))
	nelz = int(input("Number of elements in z direction, nelz = "))
	grid_type = input("XYZ grid spacing ('e' for equdistant or 'c' for "
			"Chebyshev, e.g. 'eec'): ")
	grid_function = [{'e': equidistant_grid, 'c': chebyshev_grid}[i] for i
			in grid_type]
	ind_order = input("Element index order ('xyz', 'xzy', 'yxz', 'yzx',"
			" 'zxy' or 'zyx): ")
	num_ind_order = [{'x': 0, 'y': 1, 'z': 2}[i] for i in ind_order]

	print("\nPlease specify boundary conditions ('n' for no-slip,"
			" 'f' for free-slip and 'p' for periodic)")
	bc_text = ("{0} ({1}-direction) boundary condition: ")
	bc = [input(bc_text.format("West", "x")), input(bc_text.format("East",
			"x")), input(bc_text.format("South", "y")),
			input(bc_text.format("North", "y")),
			input(bc_text.format("Bottom", "z")),
			input(bc_text.format("Top", "z"))]
	bc_num = [{'f': -2, 'n': -1, 'p': 0}[i] for i in bc]
	filename = input("\nFilename: ")

	return (np.array([xmin, ymin, zmin]), np.array([xmax, ymax, zmax]),
			np.array([nelx, nely, nelz]), grid_function,
			num_ind_order, bc_num, filename)

def main():
	n = 2
	xmin, xmax, nel, grid_func, ind_order, bc, filename = read_user_input()
	f = h5py.File(filename, "w")
	write_parameters(f, n, nel[0] * nel[1]* nel[2])
	write_elements(f, n, nel, grid_func, ind_order, xmin, xmax, bc)
	f.close()

if __name__ == "__main__":
	main()
