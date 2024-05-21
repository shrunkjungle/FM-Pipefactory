"""

def test_pipe_factory():

    import pipefactory as pf
    
    Straight0 = {
        'Name': 'Straight0',
        'length': 10.0,
        'type': 'Straight',
        }

    Bend0 = {
                'Name': 'Bend0',
                'type': 'Bend',
                'param': {
                    'radius': 20.0,
                    'angle': 30.0,
                    'axis': "up_down"
                }
            }

    Straight1 = {
                'Name': 'Straight1',
                'length': 10.0,
                'type': 'Straight',
            }

    section_list = [Straight0, Bend0, Straight1]

    mesh = pf.Pipe(70.0, 3.0, section_list, ("quad", True), 5.0, 3)


def test_section_ends_rotate_down():

    import pipefactory as pf
    import numpy as np

    Straight0 = {
        'Name': 'Straight0',
        'length': 300.0,
        'type': 'Straight',
    }

    Bend2 = {
        'Name': 'Bend1',
        'type': 'Bend',
        'param': {
            'axis' : "up_down",
            'radius': 200.0,
            'angle': -90.0
        }
    }

    section_list = [Straight0, Bend2]

    mesh = pf.Pipe(70.0, 3.0, section_list, ("hex", False), 6., 3)

    np.testing.assert_almost_equal(mesh.triads[2][0], np.array([0., 0., -1.]))
    np.testing.assert_almost_equal(mesh.triads[2][1], np.array([0., 1., 0.]))
    np.testing.assert_almost_equal(mesh.triads[2][2], np.array([1., 0., 0.]))

    np.testing.assert_almost_equal(mesh.x_sec_ends[0], np.array([0., 0., 0.]))
    np.testing.assert_almost_equal(mesh.x_sec_ends[1], np.array([300., 0., 0.]))
    np.testing.assert_almost_equal(mesh.x_sec_ends[2], np.array([500., 0., -200.]))


def test_section_ends_rotate_up():

    import pipefactory as pf
    import numpy as np

    Straight0 = {
        'Name': 'Straight0',
        'length': 350.0,
        'type': 'Straight',
    }

    Bend2 = {
        'Name': 'Bend1',
        'type': 'Bend',
        'param': {
            'axis' : "up_down",
            'radius': 150.0,
            'angle': 90.0
        }
    }

    section_list = [Straight0, Bend2]

    mesh = pf.Pipe(70.0, 3.0, section_list, ("hex", False), 6., 3)

    np.testing.assert_almost_equal(mesh.triads[2][0], np.array([0., 0., 1.]))
    np.testing.assert_almost_equal(mesh.triads[2][1], np.array([0., 1., 0.]))
    np.testing.assert_almost_equal(mesh.triads[2][2], np.array([-1., 0., 0.]))

    np.testing.assert_almost_equal(mesh.x_sec_ends[0], np.array([0., 0., 0.]))
    np.testing.assert_almost_equal(mesh.x_sec_ends[1], np.array([350., 0., 0.]))
    np.testing.assert_almost_equal(mesh.x_sec_ends[2], np.array([500., 0., 150.]))


def test_section_ends_rotate_right():

    import pipefactory as pf
    import numpy as np

    Straight0 = {
        'Name': 'Straight0',
        'length': 350.0,
        'type': 'Straight',
    }

    Bend1 = {
        'Name': 'Bend1',
        'type': 'Bend',
        'param': {
            'axis' : "left_right",
            'radius': 200.0,
            'angle': 90.0
        }
    }

    section_list = [Straight0, Bend1]

    mesh = pf.Pipe(70.0, 3.0, section_list, ("hex", False), 6., 3)

    np.testing.assert_almost_equal(mesh.triads[2][0], np.array([0., 1., 0.]))
    np.testing.assert_almost_equal(mesh.triads[2][1], np.array([-1., 0., 0.]))
    np.testing.assert_almost_equal(mesh.triads[2][2], np.array([0., 0., 1.]))

    np.testing.assert_almost_equal(mesh.x_sec_ends[0], np.array([0., 0., 0.]))
    np.testing.assert_almost_equal(mesh.x_sec_ends[1], np.array([350., 0., 0.]))
    np.testing.assert_almost_equal(mesh.x_sec_ends[2], np.array([550., 200., 0.]))

def test_section_ends_rotate_left():

    import pipefactory as pf
    import numpy as np

    Straight0 = {
        'Name': 'Straight0',
        'length': 350.0,
        'type': 'Straight',
    }

    Bend1 = {
        'Name': 'Bend1',
        'type': 'Bend',
        'param': {
            'axis' : "left_right",
            'radius': 200.0,
            'angle': -90.0
        }
    }

    section_list = [Straight0, Bend1]

    mesh = pf.Pipe(70.0, 3.0, section_list, ("hex", False), 6., 3)

    np.testing.assert_almost_equal(mesh.triads[2][0], np.array([0., -1., 0.]))
    np.testing.assert_almost_equal(mesh.triads[2][1], np.array([1., 0., 0.]))
    np.testing.assert_almost_equal(mesh.triads[2][2], np.array([0., 0., 1.]))

    np.testing.assert_almost_equal(mesh.x_sec_ends[0], np.array([0., 0., 0.]))
    np.testing.assert_almost_equal(mesh.x_sec_ends[1], np.array([350., 0., 0.]))
    np.testing.assert_almost_equal(mesh.x_sec_ends[2], np.array([550., -200., 0.]))

    """