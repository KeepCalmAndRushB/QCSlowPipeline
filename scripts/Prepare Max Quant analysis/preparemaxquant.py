import xml.etree.ElementTree as ET


def clear_parameter_tree(xml_tree):
    """Returns an xml tree object from MaxQuant's mqpar.xml that has some parts missing."""
    root = xml_tree.getroot()
    for name in ['filePaths', 'experiments', 'fractions', 'ptms', 'paramGroupIndices', 'commonChannel']:
        Child = next(root.iter(name))
        for _ in range(2):
            for short in Child:
                Child.remove(short)


def fill_parameter_tree(tree, filename):
    """Modifies an xml tree object from MaxQuant's mqpar.xml file to include a raw filename to process."""
    data = [
        ('filePaths', 'string', filename),
        ('experiments', 'short', filename),
        ('fractions', 'short', '32767'),
        ('ptms', 'boolean', 'False'),
        ('paramGroupIndices', 'int', '0'),
        ('commonChannel', 'string', '')
    ]

    for child, tag, text in data:
        element = ET.Element(tag)
        element.text = text
        next(tree.iter(child)).append(element)

def main(rawfilename, template_filename, num_threads, parameterfilename):

    tree = ET.parse(template_filename)
    clear_parameter_tree(tree)

    fill_parameter_tree(tree, filename=rawfilename)

    numThreads = next(tree.iter('numThreads'))
    numThreads.text = str(num_threads)
    # ET.dump(numThreads)

    tree.write(parameterfilename)


if __name__ == '__main__':
    main(rawfilename, template_filename, num_threads, parameterfilename)