from pymol import cmd, plugins
from pymol.plugins import addmenuitemqt
import tkinter as tk
from tkinter import simpledialog

def list_contacts(selection1, selection2, cutoff=3.5, mode='all'):
    """
    Lists contacts between two selections within a specified distance cutoff.
    """
    def callback(atom1, atom2, distance, mode=mode):
        print(f"Contact found: {atom1.index} {atom1.name} of {atom1.resn}{atom1.resi} to {atom2.index} {atom2.name} of {atom2.resn}{atom2.resi}, distance {distance:.2f} Å")

    cmd.find_pairs(f"({selection1}) within {cutoff} of ({selection2})", f"({selection2}) within {cutoff} of ({selection1})", callback=callback, mode=mode)

def plugin_app():
    parent = tk.Tk()
    parent.withdraw()
    selection1 = simpledialog.askstring("Input", "Enter first selection:", parent=parent)
    selection2 = simpledialog.askstring("Input", "Enter second selection:", parent=parent)
    cutoff = simpledialog.askfloat("Input", "Enter distance cutoff:", parent=parent)

    if selection1 and selection2 and cutoff is not None:
        list_contacts(selection1, selection2, cutoff=cutoff)
    parent.destroy()

def __init_plugin__(app=None):
    addmenuitemqt('List Contacts', plugin_app)
