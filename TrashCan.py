X = make_mesh(xmin, xmax, dx)
    for i in range(ngc):
        X.pop(0)
        X.pop(-1)

plt.plot(X, u, marker='.')