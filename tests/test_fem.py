def test_B_matrix():

    import numpy as np
    import pipefactory as pf

    le = pf.LinearElasticity(reduced = True)

    B = le.computeBMatrix(le.gauss.points[0].dN)

    ## Case 1 - Uniform Strain in X
    u = np.zeros((24,))
    id = le.index_u[[1,2,5,6]]
    u[id] = 1.0
    np.testing.assert_almost_equal(np.matmul(B, u), np.array([0.5, 0.0, 0.0, 0.0, 0.0, 0.0]))

    # Case 2 - Uniform Strain in Z
    u = np.zeros((24,))
    id = le.index_w[[4,5,6,7]]
    u[id] = 1.0
    np.testing.assert_almost_equal(np.matmul(B, u), np.array([0.0, 0.0, 0.5, 0.0, 0.0, 0.0]))

    # Case 3 - Uniform Strain in Y
    u = np.zeros((24,))
    id = le.index_v[[2,3,6,7]]
    u[id] = 1.0
    np.testing.assert_almost_equal(np.matmul(B, u), np.array([0.0, 0.5, 0.0, 0.0, 0.0, 0.0]))


def test_isoparametric():

    import numpy as np
    import pipefactory as pf

    mesh = pf.UnitBrick()

    le = pf.LinearElasticity(reduced = True)

    def getX(e):
        ln = e.list_of_nodes
        X = []
        for n in ln:
            X.append(mesh.nodes[n].coords)
        return np.array(X).transpose()

    invJ, detJ = le.computeInvJacobian(getX(mesh.elements[0]), le.gauss.points[0])

    np.testing.assert_almost_equal(0.125, detJ)
    np.testing.assert_almost_equal(6.0, invJ.trace())

    V = 0.0

    for i, ip in enumerate(le.gauss.points):
        _, detJ = le.computeInvJacobian(getX(mesh.elements[0]), ip)
        V += detJ * ip.wgt

    np.testing.assert_almost_equal(1.0, V)


def test_stiffness_matrix():

    import random
    import numpy as np
    import pipefactory as pf

    Straight0 = {
        'Name': 'Straight0',
        'length': 350.0,
        'type': 'Straight',
    }

    mesh = pf.Pipe(radius = 70., thickness= 3.0, section_list=[Straight0], elem_type=("hex", False), element_size=6.0, elements_through_thickness=3)

    le = pf.LinearElasticity(reduced = False)

    def getX(e):
        ln = e.list_of_nodes
        X = []
        for n in ln:
            X.append(mesh.nodes[n].coords)
        return np.array(X).transpose()

    for i in range(10):

        e = random.choice(mesh.elements)
        Ke = le.computeStiffnessMatrix(getX(e), E = 1.0, nu = 0.2)
        lam, v = np.linalg.eigh(Ke)

        np.testing.assert_almost_equal(len(np.where(lam < 1.e-10)[0]), 6) # Check 6 zero energy modes
        np.testing.assert_almost_equal(len(np.where(lam < -1.e-10)[0]), 0) # Check Positive Definate


