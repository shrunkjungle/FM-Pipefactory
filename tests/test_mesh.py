def test_element():

    from ..pipefactory import Element

    ele = Element([0, 1, 2], "tri", 0, 0)

    ele.change_node(2, 3)
    assert ele.list_of_nodes[2] == 3

    id, lon = ele.for_xml()
    assert id == 0
    assert lon[1] == 1