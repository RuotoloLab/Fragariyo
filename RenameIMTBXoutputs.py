"""
Take IMTBX/Grppr outputs that are individual .isotopes files in a folder named after their raw
data file and rename them to (raw name).isotopes and place all in output folder.
#author: DP
#date: 8/13/18
"""

import os
from PyQt5 import QtWidgets
import sys
CONFIG_FILE = 'config.txt'  # config file for saving last directory for fancy filedialog


def get_data(config_file):
    """
    Load folders of data using custom FileDialog class
    :param config_file: path to the config file with the initial directory for the file chooser
    :return: list of strings of full system folder paths to the folders chosen, updated input_dir
    """
    input_dir = get_last_dir(config_file)

    app = QtWidgets.QApplication(sys.argv)
    ex = FileDialog(input_dir)
    ex.show()
    app.exec_()
    files = ex.selectedFiles()

    new_base_dir = os.path.dirname(files[0])
    save_config(config_file, new_base_dir)
    return files


def get_last_dir(config_file):
    """
    parse the config file for the last directory used, to use as the initial directory when
    opening the file chooser.
    :param config_file: text file with a single directory (full system path) and nothing else
    :return: (string) directory path
    """
    with open(config_file, 'r') as config:
        line = config.readline()
        return line


def save_config(config_file, new_base_dir):
    """
    Update the config file with a new directory name
    :param config_file: file path to update
    :param new_base_dir: information to save in the config file
    :return: void
    """
    with open(config_file, 'w') as config:
        config.write(new_base_dir)


class FileDialog(QtWidgets.QFileDialog):
    """
    File chooser for raw data, created after extensive searching on stack overflow
    """
    def __init__(self, input_dir, *args):
        QtWidgets.QFileDialog.__init__(self, *args)
        self.setOption(self.DontUseNativeDialog, True)
        self.setFileMode(self.DirectoryOnly)
        self.setDirectory(input_dir)

        self.tree = self.findChild(QtWidgets.QTreeView)
        self.tree.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)


# The code above allows the Renamer to be independent from DmitryPeakToolWrapper
def main():
    """
    Choose directories containing .isotopes files to rename and merge
    :return: void
    """
    selected_dirs = get_data(CONFIG_FILE)
    main_dir = os.path.split(selected_dirs[0])[0]

    for folder in selected_dirs:
        os.chdir(folder)
        files = [x for x in os.listdir(folder) if x.endswith('.isotopes') and x.startswith('func')]
        new_filename = os.path.join(main_dir, os.path.basename(folder).rstrip('.raw') + '.isotopes')

        # Rename outputs from IMAnnotatorTD_v2 (from batch_IMTBX_process) to be the raw filename rather than func001.hits
        # files = [x for x in os.listdir(folder) if x.endswith('.hits') and x.startswith('func')]
        # new_filename = os.path.join(main_dir, os.path.basename(folder).rstrip('.raw') + '.hits')

        try:
            os.rename(files[0], new_filename)
        except FileExistsError:
            print('***overwriting old file***')
            os.remove(new_filename)
            os.rename(files[0], new_filename)


if __name__ == '__main__':
    main()
