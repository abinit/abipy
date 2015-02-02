from __future__ import print_function, division, unicode_literals

from fireworks import Workflow


def parse_workflow(fws, links_dict):
    new_list = []
    for fw in fws:
        if isinstance(fw, Workflow):
            new_list.extend(fw.fws)
        else:
            new_list.append(fw)

    new_links_dict = {}
    for parent, children in links_dict.items():
        if isinstance(parent, Workflow):
            new_links_dict.update(parent.links)
            for leaf_fw_id in parent.leaf_fw_ids:
                new_links_dict[leaf_fw_id] = children
        else:
            new_links_dict[parent] = children

    # dict since the collection will be updated
    for parent, children in dict(new_links_dict).items():
        final_childrens = []
        for child in children:
            if isinstance(child, Workflow):
                new_links_dict.update(child.links)
                final_childrens.extend(child.root_fw_ids)
            else:
                final_childrens.append(child)
        new_links_dict[parent] = final_childrens

    return new_list, new_links_dict