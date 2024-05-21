import pytest
from pipefactory import find_section  # Replace 'your_module' with the actual name of your module

@pytest.mark.parametrize("section_ends, x, expected", [
    ([1, 2, 3], 0.5, 0),   # x before first section end
    ([1, 2, 3], 1.0, 0),   # x at first section end
    ([1, 2, 3], 1.5, 1),   # x in middle section
    ([1, 2, 3], 3.0, 2),   # x at last section end
    ([1, 2, 3], 3.1, None) # x beyond last section end
    # Add more cases to cover different scenarios
])
def test_find_section(section_ends, x, expected):
    assert find_section(section_ends, x) == expected

