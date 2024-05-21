import importlib.util
import pytest


def is_extra_installed(package_name):
    return importlib.util.find_spec(package_name) is not None


PARALLEL_ENABLED = is_extra_installed("mpi4py")


def test_install():
    import pipefactory as pf
    print("This test will always run.")



@pytest.mark.skipif(PARALLEL_ENABLED, reason="requires parallel feature to be disabled")
def test_serial_feature():
    print("This test will only run if the parallel feature is NOT enabled.")


@pytest.mark.skipif(not PARALLEL_ENABLED, reason="requires parallel feature to be enabled")
def test_parallel_feature():
    print("This test will only run if the parallel feature is enabled.")
