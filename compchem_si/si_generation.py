"""si_generation
Module to take Gaussian output files and generate SI pdf pages

"""
import os

# pylint:disable=import-error
import subprocess

import cclib
import numpy as np

# from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
)

no_to_symb = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Nh",
    114: "Fl",
    115: "Mc",
    116: "Lv",
    117: "Ts",
    118: "Og",
}


def _get_geom(data):
    """Gets the geometry of the optimized structure from cclib parsed data

    Args:
        data: cclib data

    Returns:
        geometry as a list of lists with symbols and xyz coordinates

    """
    opt_geom = data.atomcoords[data.OPT_DONE - 1]
    atom_symb = [no_to_symb[no] for no in data.atomnos]
    geom_array = np.c_[np.array(atom_symb), opt_geom]
    geom_list = list(list(el) for el in geom_array)
    return geom_list


def _write_image(data, log_file):
    """Create an .xyz file for the final geometry and use that and openbabel to create an image

    Args:
        data: data of log_file parsed by cclib.io.ccread
        loog_file: name of log file

    Returns:
        None, but writes image to file

    """
    if not os.path.exists("png_files/"):
        os.mkdir("png_files")
    if not os.path.exists("xyz_files/"):
        os.mkdir("xyz_files")

    xyz_file = log_file.replace(".log", ".xyz")
    png_file = log_file.replace(".log", ".png")
    cclib.io.ccwrite(data, "xyz", f"xyz_files/{xyz_file}", indices=data.OPT_DONE - 1)
    subprocess.check_output(
        [
            "obabel",
            f"xyz_files/{xyz_file}",
            "-O",
            f"png_files/{png_file}",
            "-xp",
            "600",
            "-x0",
            "molfile",
        ],
        text=True,
    )


def parse_log_file(log_file):
    """Take a log file and get the geometry and energy values

    Args:
        log_file: the name of a .log file

    Returns:
        the energy and geometry of the molecule in the .log file, and an image written to log_file.png

    """
    data = cclib.io.ccread(log_file)
    geom_array = _get_geom(data)
    energy = data.scfenergies[data.OPT_DONE - 1]
    _write_image(data, log_file)
    return energy, geom_array


def constuct_si(log_file_list="dir", out_name="SupplementaryInformation.pdf"):
    """Construct the Supplementary Information

    Args:
        log_file_list: either 'dir' to run for all .log files in the current directory or a list of .log files to run
        out_name: name for the resulting document including the .pdf extension

    Returns:
        None, but writes generated SI to out_name
    """
    if log_file_list == "dir":
        log_file_list = [file for file in os.listdir() if file.endswith(".log")]
    doc = SimpleDocTemplate(out_name)
    story = []
    for log_file in log_file_list:
        story = construct_si_page(log_file, story)
    doc.build(story)


def construct_si_page(log_file, story):
    """For a given molecule, construct the SI page

    Args:
        log_file: log file name
        story: list containing platypus flowables from reportlab for the document being built

    Returns:
        the story input updated with the current page being built

    """
    energy, geom = parse_log_file(log_file)
    story = _supporting_info_page(story, energy, geom, log_file)
    return story


def _supporting_info_page(story, energy, geom, log_file):
    """Construct a SI page for a given molecule

    Args:
        story: list containing platypus flowables from reportlab for the document being built
        Energy: the energy for the molecule in the .log file
        geom: list[list] the geometry of the molecule
        log_file: the log file name

    Returns:
        the story input updated with the current page being built
    """

    p = Paragraph("<font size=20>Structure</font>")
    story.append(p)
    story.append(Spacer(1, 0.2 * inch))

    I = Image(f'png_files/{log_file.replace(".log", ".png")}')
    I.drawHeight = 2 * inch
    I.drawWidth = 2 * inch
    story.append(I)
    p = Paragraph(f"<font size=20>Energy</font>: {energy}")
    story.append(p)
    story.append(Spacer(1, 0.2 * inch))
    p = Paragraph("<font size=20>Geometry</font>")
    story.append(p)
    story.append(Spacer(1, 0.2 * inch))
    t = Table(geom)
    story.append(t)
    story.append(PageBreak())
    return story


if __name__ == "__main__":
    log_files = [file for file in os.listdir() if file.endswith(".log")]
    constuct_si(log_files)
