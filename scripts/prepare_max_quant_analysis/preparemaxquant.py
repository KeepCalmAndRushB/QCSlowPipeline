import xml.etree.ElementTree as ET
import argparse


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

def main(rawfilename, template_filename, analysis_directory, num_threads):

    tree = ET.parse(template_filename)
    clear_parameter_tree(tree)

    fill_parameter_tree(tree, filename=rawfilename)

    numThreads = next(tree.iter('numThreads'))
    numThreads.text = str(num_threads)
    # ET.dump(numThreads)

    tree.write(analysis_directory + "\mqpar.xml")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Prepare MaxQuant Analysis.")
    parser.add_argument('rawfilename', default='D:/1/test/test.raw', help="Raw-File Path in analysis directory")
    parser.add_argument('template_filename', default='C:/MQ/mqpar.xml', help="MQ parameter file complete path")
    parser.add_argument('directory', default='D:/1/test', help="The Folder where the MQ analysis will take place")
    parser.add_argument('threads', default=2, help="Number of threads used for analysis (default=2)")
    args = parser.parse_args()
    main(rawfilename=args.rawfilename, template_filename=args.template_filename, analysis_directory=args.directory, num_threads=args.threads)

