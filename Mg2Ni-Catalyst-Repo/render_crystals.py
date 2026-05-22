from ase.io import read, write
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms

# --- CONFIGURATION ---
# Change this string to rotate the crystal structures in the output images.
# Format is 'Xx,Yy,Zz' in degrees (e.g., '30x,-70y,15z').
VIEW_ROTATION = '15x,-30y,0z'
# ---------------------

import numpy as np

def draw_axes_triad(ax, rotation_str):
    R = np.eye(3)
    if rotation_str:
        for r in rotation_str.split(','):
            axis = r[-1]
            angle = float(r[:-1]) * np.pi / 180.0
            c, s = np.cos(angle), np.sin(angle)
            if axis == 'x': r_mat = np.array([[1,0,0], [0,c,s], [0,-s,c]])
            elif axis == 'y': r_mat = np.array([[c,0,-s], [0,1,0], [s,0,c]])
            elif axis == 'z': r_mat = np.array([[c,s,0], [-s,c,0], [0,0,1]])
            R = np.dot(R, r_mat)
            
    # ASE's projection aligns Z out of the screen. We project X,Y,Z vectors
    # We will place the triad in the bottom left corner dynamically
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    origin = np.array([xlim[0] + 0.1*(xlim[1]-xlim[0]), ylim[0] + 0.1*(ylim[1]-ylim[0])])
    length = 0.1 * (xlim[1] - xlim[0])
    
    # Project basis vectors
    basis = np.eye(3) * length
    proj = np.dot(basis, R)[:, :2]  # Keep only X and Y coordinates on screen
    
    colors = ['r', 'g', 'b']
    labels = ['X', 'Y', 'Z']
    for i in range(3):
        ax.annotate('', xy=origin + proj[i], xytext=origin,
                    arrowprops=dict(arrowstyle="->", color=colors[i], lw=2))
        ax.text(origin[0] + proj[i,0]*1.2, origin[1] + proj[i,1]*1.2, labels[i],
                color=colors[i], fontsize=12, fontweight='bold', ha='center', va='center')

def render(filepath, outpath, title, format_type, rotation=VIEW_ROTATION):
    atoms = read(filepath, format=format_type)
    
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_atoms(atoms, ax, radii=0.8, rotation=rotation)
    ax.set_axis_off()
    
    # Draw the coordinate axes
    draw_axes_triad(ax, rotation)
    ax.set_title(title, fontsize=16, pad=15)
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"Rendered {outpath}")
def render_all_in_one(outpath, rotation=VIEW_ROTATION):
    files = [
        ('mg2nih4-28.pwi', '(a) Pristine Mg$_2$NiH$_4$', 'espresso-in'),
        ('mg2nih4_Nb.in', '(b) Mg$_2$NiH$_4$:Nb', 'espresso-in'),
        ('mg2nih4_NbFe.in', '(c) Mg$_2$NiH$_4$:Nb,Fe', 'espresso-in')
    ]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    unique_symbols = set()
    
    for i, (ax, (filepath, title, format_type)) in enumerate(zip(axes, files)):
        atoms = read(filepath, format=format_type)
        unique_symbols.update(atoms.get_chemical_symbols())
        plot_atoms(atoms, ax, radii=0.8, rotation=rotation)
        ax.set_axis_off()
        if i == 0:
            draw_axes_triad(ax, rotation)
        ax.set_title(title, fontsize=16, pad=15)
        
    # Add legend at the bottom
    from ase.data.colors import jmol_colors
    from ase.data import atomic_numbers
    from matplotlib.lines import Line2D
    
    handles = [
        Line2D([0], [0], marker='o', color='w', 
               markerfacecolor=jmol_colors[atomic_numbers[sym]], 
               markeredgecolor='k', markersize=12, label=sym) 
        for sym in sorted(unique_symbols)
    ]
    fig.legend(handles=handles, loc='lower center', ncol=len(handles), 
               fontsize=16, frameon=False, bbox_to_anchor=(0.5, 0.20))
        
    plt.tight_layout()
    # Adjust bottom to make room for legend
    plt.subplots_adjust(bottom=0.25)
    plt.savefig(outpath, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"Rendered {outpath}")


if __name__ == '__main__':
    # System 0: Pristine Mg2NiH4
    render('mg2nih4-28.pwi', 'tex/struct_sys0.png', '(a) Pristine Mg$_2$NiH$_4$', 'espresso-in')
    # System 1: Doped Mg2NiH4 (Nb)
    render('mg2nih4_Nb.in', 'tex/struct_sys1.png', '(b) Mg$_2$NiH$_4$:Nb', 'espresso-in')

    # System 2: Doped Mg2NiH4 (Nb, Fe) - using .in because fast_sweep didn't finish an .out file completely for this one yet
    render('mg2nih4_NbFe.in', 'tex/struct_sys2.png', '(c) Mg$_2$NiH$_4$:Nb,Fe', 'espresso-in')

    # Combined figure
    render_all_in_one('tex/struct_combined.png')
